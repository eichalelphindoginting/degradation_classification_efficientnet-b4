// 1. Konfigurasi Penelitian
var CONFIG = {
  firstYear: 2024,
  firstMonth: 1,
  lastYear: 2026,
  lastMonth: 7,
  crs: 'EPSG:32748',
  scaleM: 10,
  patchSizePx: 128,
  cloudScoreThreshold: 0.60,
  topoIlluminationThreshold: 0.15,
  ndviDegradedThreshold: 0.45,
  ndviHealthyThreshold: 0.65,
  minDegFraction: 0.02,
  maxDegHealthy: 0.005,
  minHealthyFraction: 0.50,
  minClearFraction: 0.80,
  minForestFraction: 0.50,
  minValidForestFraction: 0.80,
  forestPrefilterScaleM: 30,
  forestPrefilterMargin: 0.05,
  tileScale: 16,
  metadataExportMode: 'yearly',
  exportMetadata: true,
  exportImages: true,
  showPreview: false,
  imageNoData: -9999,
  exportFolder: 'TNBBS_CNN_READY_DATASET'
};

//Ukuran Grid dan Kolom Metadata yang akan dibuat
var GRID_STEP_M = CONFIG.patchSizePx * CONFIG.scaleM;
var UTM_PROJECTION = ee.Projection(CONFIG.crs);
var EXPORT_TRANSFORM = [CONFIG.scaleM, 0, 0, 0, -CONFIG.scaleM, 10000000];
var METADATA_SELECTORS = [
  'grid_id', 'target_month', 'year', 'month',
  'label', 'label_name',
  'degraded_fraction', 'healthy_fraction',
  'clear_fraction', 'forest_fraction', 'valid_forest_fraction',
  'mean_ndvi', 'hansen_fraction', 'hansen_available',
  'lon', 'lat', 'x_utm', 'y_utm'
];


// 2. Batas Wilayah Penelitian Menggunakan FAO GAUL
// Taman Nasional Bukit Barisan Selatan (TNBBS)
var tnbbs = ee.FeatureCollection('WCMC/WDPA/current/polygons')
  .filter(ee.Filter.or(
    ee.Filter.stringContains('NAME', 'Bukit Barisan Selatan'),
    ee.Filter.stringContains('ORIG_NAME', 'Bukit Barisan Selatan')
  ));

// Provinsi Lampung
var lampung = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filter(ee.Filter.eq('ADM1_NAME', 'Lampung'));
var aoi = tnbbs.geometry().intersection(lampung.geometry(), 100).simplify(100);
var exportRegion = aoi.bounds(100).buffer(GRID_STEP_M);


// 3. Forest Mask Untuk Pelabelan dan Statistik ( Tidak Diterapkan Pada Input CNN )
var hansen = ee.Image('UMD/hansen/global_forest_change_2025_v1_13');
var dem = ee.Image('USGS/SRTMGL1_003').select('elevation');

// Prefilter untuk menentukan area hutan, perairan, dan ketinggian di atas pantai (20m)
// Daerah Hutan >= 30% Canopy cover + Tidak ada perairan + Ketinggian > 20m 
var hansenLand = hansen.select('treecover2000').gte(30)
  .and(hansen.select('datamask').eq(1));
var notWater = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
  .select('occurrence').gte(50).unmask(0).not();
var aboveCoast = dem.gte(20);

var forestMask = hansenLand.and(notWater).and(aboveCoast)
  .rename('forest_mask');


// 4. Pembentukan Grid 1280 x 1280 m serta membentuk grid_id unik pada setiap bulan
var rawGrid = aoi.coveringGrid({
  proj: CONFIG.crs,
  scale: GRID_STEP_M
});

var grid = rawGrid.map(function(feature) {
  var centroidUtm = feature.geometry().centroid(1, UTM_PROJECTION);
  var coordinates = centroidUtm.coordinates();
  var xUtm = ee.Number(coordinates.get(0)).round();
  var yUtm = ee.Number(coordinates.get(1)).round();
  var gridId = ee.String('grid_').cat(xUtm.format('%.0f'))
    .cat('_').cat(yUtm.format('%.0f'));

  return feature.set({
    grid_id: gridId,
    x_utm: xUtm,
    y_utm: yUtm
  });
});

// Prefilter Grid Hutan 30m
var gridForReduction = grid.limit(100000);
// Menghitung nilai fraksi hutan pada setiap grid
var forestPrefilterImage = forestMask.unmask({
  value: 0,
  sameFootprint: false
}).rename('forest_prefilter_fraction').clip(exportRegion);

var rawForestPrefilter = forestPrefilterImage.reduceRegions({
  collection: gridForReduction,
  reducer: ee.Reducer.mean().setOutputs(['forest_prefilter_fraction']),
  scale: CONFIG.forestPrefilterScaleM,
  crs: CONFIG.crs,
  tileScale: CONFIG.tileScale
});

var prefilteredGrid = rawForestPrefilter
  .filter(ee.Filter.notNull(['forest_prefilter_fraction']))
  .filter(ee.Filter.gte(
    'forest_prefilter_fraction',
    CONFIG.minForestFraction - CONFIG.forestPrefilterMargin
  ))
  .limit(100000);

// Fraksi Final Hutan yang dihitung menggunakan resolusi 10 m
var forestFractionImage = forestMask.unmask({
  value: 0,
  sameFootprint: false
}).rename('forest_fraction').clip(exportRegion);

var rawGridWithForest10m = forestFractionImage.reduceRegions({
  collection: prefilteredGrid,
  reducer: ee.Reducer.mean().setOutputs(['forest_fraction']),
  scale: CONFIG.scaleM,
  crs: CONFIG.crs,
  tileScale: CONFIG.tileScale
});

var analysisGrid = rawGridWithForest10m
  .filter(ee.Filter.notNull(['forest_fraction']))
  .filter(ee.Filter.gte('forest_fraction', CONFIG.minForestFraction))
  .limit(100000);


// 5. Masking untuk Dataset Citra Sentinel-2
var terrain = ee.Algorithms.Terrain(dem);
var slopeRad = terrain.select('slope').multiply(Math.PI / 180);
var aspectRad = terrain.select('aspect').multiply(Math.PI / 180);

// Fungsi Filtrasi SCL (Cloud Mask) untuk membuang piksel yang 
// tertutup awan, bayangan, salju, dan wilayah yang tidak memenuhi syarat
function addQualityMask(image) {
  var scl = image.select('SCL');
  var badScl = scl.eq(1).or(scl.eq(2)).or(scl.eq(3)).or(scl.eq(6))
    .or(scl.eq(8)).or(scl.eq(9)).or(scl.eq(10)).or(scl.eq(11));
  var cloudScoreMask = image.select('cs_cdf')
    .gte(CONFIG.cloudScoreThreshold);

// Fungsi Filtrasi Topografi untuk mengetahu efek pencahayaan matahari pada permukaan bumi
// Digunakan untuk menghitung nilai NDVI yang representatif
  var solarZenithRad = ee.Image.constant(
    ee.Number(image.get('MEAN_SOLAR_ZENITH_ANGLE')).multiply(Math.PI / 180)
  );
  var solarAzimuthRad = ee.Image.constant(
    ee.Number(image.get('MEAN_SOLAR_AZIMUTH_ANGLE')).multiply(Math.PI / 180)
  );

  var cosineIncidence = solarZenithRad.cos().multiply(slopeRad.cos())
    .add(solarZenithRad.sin().multiply(slopeRad.sin())
      .multiply(solarAzimuthRad.subtract(aspectRad).cos()));
  var topographicMask = cosineIncidence.gt(
    CONFIG.topoIlluminationThreshold
  );

  return image
    .updateMask(cloudScoreMask.and(badScl.not()).and(topographicMask))
    .select(['B2', 'B3', 'B4', 'B8'])
    .copyProperties(image, ['system:time_start']);
}

// Mengambil Nilai Cloud Score+ untuk lokasi dan rentang tanggal yang sama 
function getMonthlyComposite(startDate, endDate) {
  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(exportRegion)
    .filterDate(startDate, endDate);
  var cloudScorePlus = ee.ImageCollection(
    'GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED'
  )
    .filterBounds(exportRegion)
    .filterDate(startDate, endDate);

  return s2.linkCollection(cloudScorePlus, ['cs_cdf'])
    .map(addQualityMask)
    .median();
}

// Penggunaan Hansen untuk mengetahui data loss, penggunaan hansen maksimal pada tahun 2025
// karena data loss terbaru dari Hansen berakhir pada tahun 2025
function getHansenMetadataImage(targetYear) {
  var lastLossCode = Math.min(Math.max(targetYear - 2000, 0), 25);
  return hansen.select('lossyear').gte(1)
    .and(hansen.select('lossyear').lte(lastLossCode))
    .and(forestMask)
    .rename('hansen_raw')
    .unmask({value: 0, sameFootprint: false});
}


// 6. Pemrosesan Penggunaan Mask, NDVI, dan Metadata Untuk Setiap Bulan
function processMonth(year, month) {
  var startDate = ee.Date.fromYMD(year, month, 1);
  var endDate = startDate.advance(1, 'month');
  var targetMonth = year + '_' + pad2(month);

  // Komposit bulanan tidak dibatasi oleh AOI namun untuk pelabelan 
  // dibatasi oleh AOI atau wilayah yang dipakai penelitian
  var monthlySource = getMonthlyComposite(startDate, endDate);
  var monthlyForStats = monthlySource.clip(aoi);

  // Perhitungan NDVI 
  var ndvi = monthlyForStats.normalizedDifference(['B8', 'B4'])
    .rename('NDVI');

  // Piksel yang bersih dari awan dan tidak tertutup topografi
  var clearMask = ndvi.mask().rename('clear_raw').unmask({
    value: 0,
    sameFootprint: false
  });
  // Piksel yang berada dalam wilayah hutan
  var validForestMask = ndvi.mask().and(forestMask);
  // Logika penentuan pelabelan piksel berdasarkan NDVI
  // Terdegradasi NDVI < 0.45
  var degradedMask = ndvi.lt(CONFIG.ndviDegradedThreshold)
    .and(validForestMask)
    .rename('deg_raw')
    .unmask({value: 0, sameFootprint: false});

  // Sehat NDVI >= 0.65
  var healthyMask = ndvi.gte(CONFIG.ndviHealthyThreshold)
    .and(validForestMask)
    .rename('healthy_raw')
    .unmask({value: 0, sameFootprint: false});
  
  // Hutan valid
  var validForestRaw = validForestMask.rename('valid_forest_raw').unmask({
    value: 0,
    sameFootprint: false
  });
  // Nilai rata-rata NDVI
  var meanNdvi = ndvi.updateMask(validForestMask).rename('mean_ndvi');
  // Data loss dari Hansen
  var hansenRaw = getHansenMetadataImage(year);

  // Menggabungkan 6 band statistik menjadi 1 gambar
  var statsImage = ee.Image.cat([
    degradedMask,
    healthyMask,
    validForestRaw,
    clearMask,
    hansenRaw,
    meanNdvi
  ]);

  // Fungsi untuk menghitung nilai rata-rata dari setiap band pada setiap grid
  var rawStats = statsImage.reduceRegions({
    collection: analysisGrid,
    reducer: ee.Reducer.mean(),
    scale: CONFIG.scaleM,
    crs: CONFIG.crs,
    tileScale: CONFIG.tileScale
  }).filter(ee.Filter.notNull([
    'deg_raw', 'healthy_raw', 'valid_forest_raw',
    'clear_raw', 'hansen_raw', 'mean_ndvi', 'forest_fraction'
  ]));

  // Fungsi untuk menghitung fraksi dan memfilter grid berdasarkan proporsi piksel
  var patchStats = rawStats.map(function(feature) {
    var degRaw = ee.Number(feature.get('deg_raw'));
    var healthyRaw = ee.Number(feature.get('healthy_raw'));
    var validForest = ee.Number(feature.get('valid_forest_raw'));
    var forestFraction = ee.Number(feature.get('forest_fraction'));
    var safeValidForest = validForest.max(0.000001);
    var safeForest = forestFraction.max(0.000001);
  
  // Rumus untuk menghitung fraksi dan menghindari pembagian dengan nol
  // Mengembalikan 1 feature untuk setiap grid
    return feature.set({
      degraded_fraction: degRaw.divide(safeValidForest),
      healthy_fraction: healthyRaw.divide(safeValidForest),
      clear_fraction: ee.Number(feature.get('clear_raw')),
      valid_forest_fraction: validForest.divide(safeForest),
      hansen_fraction: ee.Number(feature.get('hansen_raw')).divide(safeForest)
    });
  });

  // Fungsi untuk menentukan patch yang valid berdasarkan 
  var validPatches = patchStats
    .filter(ee.Filter.gte('clear_fraction', CONFIG.minClearFraction))
    .filter(ee.Filter.gte('forest_fraction', CONFIG.minForestFraction))
    .filter(ee.Filter.gte(
      'valid_forest_fraction', CONFIG.minValidForestFraction
    ));
  
  // Fungsi untuk menentukan patch yang terdegradasi
  var degradedPatches = validPatches
    .filter(ee.Filter.gte('degraded_fraction', CONFIG.minDegFraction))
    .map(function(feature) {
      return feature.set({label: 1, label_name: 'degradasi'});
    });
    
  // Fungsi untuk menentukan patch yang sehat
  var healthyPatches = validPatches
    .filter(ee.Filter.lte('degraded_fraction', CONFIG.maxDegHealthy))
    .filter(ee.Filter.gte('healthy_fraction', CONFIG.minHealthyFraction))
    .map(function(feature) {
      return feature.set({label: 0, label_name: 'hutan_sehat'});
    });
  
  // Menggabungkan patch yang terdegradasi dan sehat
  var labeledPatches = degradedPatches.merge(healthyPatches);
  // Menentukan apakah data Hansen tersedia untuk tahun tersebut
  var hansenAvailable = year <= 2025 ? 1 : 0;

  // Mengambil metadata dari setiap patch
  var metadata = labeledPatches.map(function(feature) {
    var centroidLonLat = feature.geometry().centroid(1)
      .transform('EPSG:4326', 1);
    var lonLat = centroidLonLat.coordinates();

    // Mengembalikan feature untuk membentuk metada akhir yang akan digunakan
    return ee.Feature(null, {
      grid_id: feature.get('grid_id'),
      target_month: targetMonth,
      year: year,
      month: month,
      label: feature.get('label'),
      label_name: feature.get('label_name'),
      degraded_fraction: feature.get('degraded_fraction'),
      healthy_fraction: feature.get('healthy_fraction'),
      clear_fraction: feature.get('clear_fraction'),
      forest_fraction: feature.get('forest_fraction'),
      valid_forest_fraction: feature.get('valid_forest_fraction'),
      mean_ndvi: feature.get('mean_ndvi'),
      hansen_fraction: feature.get('hansen_fraction'),
      hansen_available: hansenAvailable,
      lon: ee.Number(lonLat.get(0)),
      lat: ee.Number(lonLat.get(1)),
      x_utm: feature.get('x_utm'),
      y_utm: feature.get('y_utm')
    });
  });

  // Fungsi untuk menyiapkan citra RGB yang akan digunakan pada proses pelatihan model
  // Citra RGB diubah menjadi format Int16 agar ukurannya lebih kecil
  // Fungsi ini hanya digunakan pada tahun 2024, 2025, dan 2026
  var rgbForTraining = monthlySource.select(['B4', 'B3', 'B2'])
    .clamp(0, 10000)
    .round()
    .toInt16()
    .unmask({value: CONFIG.imageNoData, sameFootprint: false})
    .clip(exportRegion);

  return {
    targetMonth: targetMonth,
    metadata: metadata,
    rgbImage: rgbForTraining,
    ndvi: ndvi,
    degradedPixels: ndvi.lt(CONFIG.ndviDegradedThreshold)
      .and(validForestMask).selfMask()
  };
}

// 7. Fungsi untuk menyusun daftar bulan dan membentuk task
function pad2(value) {
  return value < 10 ? '0' + value : String(value);
}

// Fungsi untuk looping task sebanyak jumlah bulan yang digunakan
function buildClientMonths(firstYear, firstMonth, lastYear, lastMonth) {
  var months = [];
  var year = firstYear;
  var month = firstMonth;

  while (year < lastYear || (year === lastYear && month <= lastMonth)) {
    months.push({year: year, month: month});
    month += 1;
    if (month === 13) {
      month = 1;
      year += 1;
    }
  }
  return months;
}

var clientMonths = buildClientMonths(
  CONFIG.firstYear,
  CONFIG.firstMonth,
  CONFIG.lastYear,
  CONFIG.lastMonth
);

var metadataByYear = {};
var allMetadata = [];
var previewResult = null;

for (var i = 0; i < clientMonths.length; i++) {
  var year = clientMonths[i].year;
  var month = clientMonths[i].month;
  var result = processMonth(year, month);

  if (!metadataByYear[year]) {
    metadataByYear[year] = [];
  }
  metadataByYear[year].push(result.metadata);
  allMetadata.push(result.metadata);

  if (CONFIG.exportMetadata && CONFIG.metadataExportMode === 'monthly') {
    Export.table.toDrive({
      collection: result.metadata,
      description: 'metadata_patches_' + result.targetMonth,
      folder: CONFIG.exportFolder,
      fileNamePrefix: 'metadata_patches_' + result.targetMonth,
      fileFormat: 'CSV',
      selectors: METADATA_SELECTORS
    });
  }

  if (CONFIG.exportImages) {
    Export.image.toDrive({
      image: result.rgbImage,
      description: 'mosaic_multiband_' + result.targetMonth,
      folder: CONFIG.exportFolder,
      fileNamePrefix: 'mosaic_multiband_' + result.targetMonth,
      region: exportRegion,
      crs: CONFIG.crs,
      crsTransform: EXPORT_TRANSFORM,
      maxPixels: 1e13,
      fileFormat: 'GeoTIFF',
      formatOptions: {
        cloudOptimized: true,
        noData: CONFIG.imageNoData
      }
    });
  }

  if (CONFIG.showPreview && i === 0) {
    previewResult = result;
  }
}


// 8. Ekspor Metadata Tahunan
if (CONFIG.exportMetadata && CONFIG.metadataExportMode === 'yearly') {
  for (var exportYear = CONFIG.firstYear;
       exportYear <= CONFIG.lastYear;
       exportYear++) {
    var yearlyMetadata = ee.FeatureCollection(
      metadataByYear[exportYear]
    ).flatten();

    Export.table.toDrive({
      collection: yearlyMetadata,
      description: 'metadata_patches_' + exportYear,
      folder: CONFIG.exportFolder,
      fileNamePrefix: 'metadata_patches_' + exportYear,
      fileFormat: 'CSV',
      selectors: METADATA_SELECTORS
    });
  }
}

if (CONFIG.exportMetadata && CONFIG.metadataExportMode === 'all') {
  var masterMetadata = ee.FeatureCollection(allMetadata).flatten();
  Export.table.toDrive({
    collection: masterMetadata,
    description: 'metadata_patches_ALL_MONTHS',
    folder: CONFIG.exportFolder,
    fileNamePrefix: 'metadata_patches_ALL_MONTHS',
    fileFormat: 'CSV',
    selectors: METADATA_SELECTORS
  });
}


// 9. Preview Opsional dan Ringkasan Client-Side
// Preview hanya bulan pertama agar Map tidak membangun banyak graph.
if (CONFIG.showPreview && previewResult !== null) {
  Map.centerObject(aoi, 9);
  Map.addLayer(
    previewResult.rgbImage.updateMask(
      previewResult.rgbImage.neq(CONFIG.imageNoData)
    ),
    {bands: ['B4', 'B3', 'B2'], min: 0, max: 3000},
    'RGB ' + previewResult.targetMonth
  );
  Map.addLayer(
    previewResult.ndvi,
    {min: 0, max: 0.9, palette: ['white', 'yellow', 'green']},
    'NDVI ' + previewResult.targetMonth
  );
  Map.addLayer(
    previewResult.degradedPixels,
    {palette: ['red']},
    'Pixel NDVI < 0.45 ' + previewResult.targetMonth
  );
}

var metadataTaskCount = CONFIG.metadataExportMode === 'monthly'
  ? clientMonths.length
  : (CONFIG.metadataExportMode === 'yearly'
    ? CONFIG.lastYear - CONFIG.firstYear + 1
    : 1);
var imageTaskCount = CONFIG.exportImages ? clientMonths.length : 0;

print('Periode: 2024_01 sampai 2026_07 (' + clientMonths.length + ' bulan)');
print('Task dibuat: ' + imageTaskCount + ' GeoTIFF RGB + ' +
  (CONFIG.exportMetadata ? metadataTaskCount : 0) + ' CSV metadata');
print('Spatial split 70/15/15 dilakukan di Python berdasarkan grid_id.');
