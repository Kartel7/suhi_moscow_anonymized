/** Algorithm for Landsat-5 & 8 Processing, atmospheric correction, summer SUHI calculation and clustering,
    some statistical calculations (mean LST, SUHI area) and export **/

/** Main used parameters:
 *
 * GEOMETRY:
 *
 * var roi: ROI (Big Moscow)
 * var rgnp: bounds for image search (Moscow center boundaries here)
 * var rgnperc: bounds for ndvi 10 and 90 percentiles calculation (Moscow boundaries here)
 * var agglomeration: parts of agglomeration inside ROI for analysis and statistics calculation (Big Moscow parts here)
 *     includes fields:
 *                      'agl': name of every agglomeration part
 * var lcz: polygons of local climate zones for statistics and clustering
 *     includes fields:
 *                      'id': unique id for each polygon
 *                      'lcz': lcz class
 *                      'name': unique name for each polygon
 *                      'rus': unique name for each polygon in the second language
 * 
 * var rural: polygons of rural areas inside ROI for SUHI calculation
 * 
 * OTHER:
 * 
 * vars L5_START_YEAR, L5_END_YEAR = ... : Last years to search for images of Landsat 5
 * vars L8_START_YEAR, L8_END_YEAR = ... : Last years to search for images of Landsat 8
 * vars example_yearL5, example_yearL8, merged_example_year = ... : Years for visualization
 * var periods: list of dictionaries with your time periods in format like [{start:y1, end:y2}, {start:y2+1, end:y3}, etc.]
**/

// Get raw data and define bounds (polygon) for image search
var rgnp = ee.FeatureCollection("projects/ee-olegkamalyuta/assets/msc").geometry();
Map.addLayer(rgnp, {color: 'red'}, 'rgnp');
Map.centerObject(rgnp, 8);

var L5_START_YEAR = 1984; // First year to search images
var L5_END_YEAR   = 2012; // Last year
var L8_START_YEAR = 2013;
var L8_END_YEAR   = 2024;

var example_yearL5 = 1999; // visualization year 1 (Landsat 5)
var example_yearL8 = 2019; // visualization year 2 (Landsat 8)

// Creating summer collection for one year
function createL5SummerCollection(year) {
  var startDate = ee.Date.fromYMD(year, 6, 1);
  var endDate   = ee.Date.fromYMD(year, 8, 31);
  
  var col = ee.ImageCollection("LANDSAT/LT05/C02/T1")
    .filterDate(startDate, endDate)
    .filterBounds(rgnp)
    .filter(ee.Filter.lte('CLOUD_COVER', 30)) // weaking the filter can lead in user memory excess (error)
    .map(function(img) { // using map to apply a function for each img in ImageCollection
       // setting "year" property for each image
       return img.set('year', year);
    });
  return col; // returning a collection of images with year property
}

// Same for Landsat 8 
function createL8SummerCollection(year) {
  var startDate = ee.Date.fromYMD(year, 6, 1);
  var endDate   = ee.Date.fromYMD(year, 8, 31);
  
  var col = ee.ImageCollection("LANDSAT/LC08/C02/T1")
    .filterDate(startDate, endDate)
    .filterBounds(rgnp)
    .filter(ee.Filter.lte('CLOUD_COVER', 30))
    .map(function(img) {
       return img.set('year', year);
    });
  return col;
}

// Forming a list of years for each satellite
var yearsListL5 = ee.List.sequence(L5_START_YEAR, L5_END_YEAR);
var yearsListL8 = ee.List.sequence(L8_START_YEAR, L8_END_YEAR);

// Creating a list from ImageCollection where each element is an one-year collection
var listOfCollectionsL5 = yearsListL5.map(function(y) {
  return createL5SummerCollection(y);
});

// listOfCollections → ee.List[ ImageCollection, ImageCollection, ... ]

// Same for Landsat 8
var listOfCollectionsL8 = yearsListL8.map(function(y) {
  return createL8SummerCollection(y);
});

// Merging collections into one using iterate
// Start from empty collection and merge with every one
var bigColL5 = ee.List(listOfCollectionsL5).iterate( // creating an iterator over the ImageCollection list for each year
  function(ic, acc) {
    // ic, acc are server objects, we convert them into ImageCollection
    var current = ee.ImageCollection(ic); // taking the output of iterator (one collection per year)
    var merged  = ee.ImageCollection(acc).merge(current); // merge it with others
    return merged;
  },
  ee.ImageCollection([])  // ininitial value
);

// Same for Landsat 8
var bigColL8 = ee.List(listOfCollectionsL8).iterate( 
  function(ic, acc) { 
    var current = ee.ImageCollection(ic);
    var merged  = ee.ImageCollection(acc).merge(current);
    return merged;
  },
  ee.ImageCollection([])
);

// bigCol is object-result of iterate, convert it to ee.ImageCollection
bigColL5 = ee.ImageCollection(bigColL5);
bigColL8 = ee.ImageCollection(bigColL8);

// Print general collection and number of images
print("Total Landsat-5 collection (summer) for all years:", bigColL5);
print("Total LS5 images number:", bigColL5.size());
print("Total Landsat-8 collection (summer) for all years:", bigColL8);
print("Total LS8 images number:", bigColL8.size());

// Count number of images for each year
function countImagesByYear(collection, yearsList) {
  return yearsList.map(function(y) {
    y = ee.Number(y);
    var subCol = collection.filterMetadata('year', 'equals', y);
    return ee.Dictionary({
      year: y,
      count: subCol.size()
    });
  });
}

// Count number of images for each year for Landsat 5
var yearCountsL5 = countImagesByYear(bigColL5, yearsListL5);

// Same for Landsat 8
var yearCountsL8 = countImagesByYear(bigColL8, yearsListL8);

// Merging results for Landsat 5 and Landsat 8
var combinedYearCounts = yearCountsL5.cat(yearCountsL8);

// Print number of images
print("Number of images by year for Landsat 5 and Landsat 8:", combinedYearCounts);


// 6) Visualize example image using example years l5 and l8
var sampleImage = bigColL5.filterMetadata('year', 'equals', example_yearL5).first()
// Map.addLayer(sampleImage, {bands: ['B3', 'B2', 'B1'], min:0, max:100}, 'Пример снимка L5');
// print(sampleImage.get('SPACECRAFT_ID'));
var sampleImage = bigColL8.filterMetadata('year', 'equals', example_yearL8).first()
// Map.addLayer(sampleImage, {bands: ['B4', 'B3', 'B2'], min:0, max:100}, 'Пример снимка L8');
// print(sampleImage.get('SPACECRAFT_ID'));

// PROCESSING BEGINNING

// Converting Digital Numbers DN to Lλ and then to Brightness Temperatures Tb for TIR bands
function dnTo_LTb_L5(image) {
  
  // Fixing constans since sometimes they're null in images' metadata
  var multPropL5 = 0.055375
  var addPropL5 = 1.18243
  var K1L5 = 607.76
  var K2L5 = 1260.56
  
    // Converting DN to Lλ for LS5 TIR band
    var L_lambda_TIR = image.select('B6')
                   .multiply(multPropL5)
                   .add(addPropL5)
                   .rename('L_lambda_TIR');
    
    // Apply formula for Tb retrieval
    var Tb_TIR = L_lambda_TIR.expression(
      'K2 / log((K1 / L) + 1)', {
        'K1': K1L5,
        'K2': K2L5,
        'L': L_lambda_TIR
      }).rename('Tb_TIR');
      image = image.addBands([L_lambda_TIR, Tb_TIR]);
      
  return image
  }
  
// Same for Landsat 8
   function dnTo_LTb_L8(image) {
     
    var multPropL8 = 0.0003342
    var addPropL8 = 0.1
    var K1L8 = 774.8853
    var K2L8 = 1321.0789
     
    // We don't use B11 since it has technical issues
    var L_lambda_TIR = image.select('B10')
                   .multiply(multPropL8)
                   .add(addPropL8)
                   .rename('L_lambda_TIR');

    var Tb_TIR = L_lambda_TIR.expression(
      'K2 / log((K1 / L) + 1)', {
        'K1': K1L8,
        'K2': K2L8,
        'L': L_lambda_TIR
      }).rename('Tb_TIR');
      image = image.addBands([L_lambda_TIR, Tb_TIR]);
  return image;
  }

  
// Applying a function to an entire collection
var radianceColL5 = bigColL5.map(dnTo_LTb_L5);
var radianceColL8 = bigColL8.map(dnTo_LTb_L8);

// Visualising Lλ and Tb (B6 for Landsat 5 and B10 for Landsat 8)
// var sampleImageRad5 = radianceColL5.filterMetadata('year', 'equals', example_yearL5).first();
// Map.addLayer(sampleImageRad5, {bands: ['L_lambda_TIR', 'Tb_TIR'], min: 0, max: 100}, 'Example L_lambda_B6 and Tb_B6');
// print('Processed collection with Lλ and Tb for LS5 B6:', radianceColL5);

// var sampleImageRad8 = radianceColL8.filterMetadata('year', 'equals', example_yearL8).first();
// Map.addLayer(sampleImageRad8, {bands: ['L_lambda_TIR', 'Tb_TIR'], min: 0, max: 100}, 'Example L_lambda_B10 and Tb_B10');
// print('Processed collection with Lλ and Tb for LS8 B10:', radianceColL8);

// Convering DN to Ltoa for optical channels
function dnToReflectanceAllBandsL5(image) {
  // Landsat 5 optical and SWIR bands
  var bands = ['B1','B2','B3','B4','B5','B7',];
  
  // Calculate Lλ for each channel using corresponding coefficients
  bands.forEach(function(band) {
    // Extracting the numeric part of the channel name (like "6" from "B6")
    var bandNum = band.slice(1);
    // Forming metadata names for multiplier and offset
    var multProp = 'REFLECTANCE_MULT_BAND_' + bandNum;
    var addProp = 'REFLECTANCE_ADD_BAND_' + bandNum;
    
    // Extracting coefficients from metadata
    var reflMult = ee.Number(image.get(multProp));
    var reflAdd = ee.Number(image.get(addProp));
    
    // Converting DN to Ltoa for chosen band
    var refl = image.select(band)
                   .multiply(reflMult)
                   .add(reflAdd)
                   .rename('L_toa_' + band);
                   
    // Correcting for solar zenith angle
    var sunElevation = ee.Number(image.get('SUN_ELEVATION'));
    var sunZenith = ee.Number(90).subtract(sunElevation);
    var sunZenithRadians = sunZenith.multiply(Math.PI).divide(180);
    var refl_sol_corr = refl.divide(sunZenithRadians.cos());
    
    // Adding received channel to an image
    image = image.addBands(refl_sol_corr);
  });
  return image;
}

// Same for Landsat 8
function dnToReflectanceAllBandsL8(image) {
  // Landsat 8 optical and SWIR bands
  var bands = ['B1','B2','B3','B4','B5', 'B6', 'B7'];
  
  bands.forEach(function(band) {
    var bandNum = band.slice(1);
    var multProp = 'REFLECTANCE_MULT_BAND_' + bandNum;
    var addProp = 'REFLECTANCE_ADD_BAND_' + bandNum;
    
    var reflMult = ee.Number(image.get(multProp));
    var reflAdd = ee.Number(image.get(addProp));
    
    var refl = image.select(band)
                   .multiply(reflMult)
                   .add(reflAdd)
                   .rename('L_toa_' + band);
                   
    var sunElevation = ee.Number(image.get('SUN_ELEVATION'));
    var sunZenith = ee.Number(90).subtract(sunElevation);
    var sunZenithRadians = sunZenith.multiply(Math.PI).divide(180);
    var refl_sol_corr = refl.divide(sunZenithRadians.cos());
    
    image = image.addBands(refl_sol_corr);
  });
  return image;
}

// Applying the function to an entire collection
var reflectanceColL5 = radianceColL5.map(dnToReflectanceAllBandsL5);
var reflectanceColL8 = radianceColL8.map(dnToReflectanceAllBandsL8);

// Visualizing VNIR+SWIR Reflectance for Landsat 5 and 8
var sampleImageReflL5 = reflectanceColL5.filterMetadata('year', 'equals', example_yearL5).first();
// Map.addLayer(sampleImageReflL5, {bands: ['L_toa_B3', 'L_toa_B2', 'L_toa_B1'], min: 0, max: 100}, 'Example L_toa L5');
// print('Processed collection with L5toa for all bands:', reflectanceColL5);

var sampleImageReflL8 = reflectanceColL8.filterMetadata('year', 'equals', example_yearL8).first();
// Map.addLayer(sampleImageReflL8, {bands: ['L_toa_B4', 'L_toa_B3', 'L_toa_B2'], min: 0, max: 100}, 'Example L_toa L8');
// print('ОProcessed collection with L8toa for all bands:', reflectanceColL8);
 
 
/*******************************
 * Cloud masking for Landsat 5
 *******************************/
function maskCloudsL5(image) {
  var qa = image.select('QA_PIXEL');

  // Bit 0: Fill - if 1 is set, there is no data
  var fillBit = 1 << 0;
  // Bit 6: Clear – if 1 is set, there is no cloud
  var clearBit = 1 << 6;

  // Checking that bit Clear = 1, and Fill = 0
  var isClear = qa.bitwiseAnd(clearBit).neq(0);  // bit 6 eq 1 => Clear
  var isFill = qa.bitwiseAnd(fillBit).neq(0);    // bit 0 eq 1 => Fill

  // Get pixels WHERE (Clear) AND (not fill)
  var mask = isClear.and(isFill.not());

  return image.updateMask(mask);
}

/*******************************
 * Cloud masking for Landsat 8
 *******************************/
function maskCloudsL8(image) {
  var qa = image.select('QA_PIXEL');

  // Bit 0: Fill
  var fillBit = 1 << 0;
  // Bit 2: Cirrus
  var cirrusBit = 1 << 2;
  // Bit 6: Clear (including Cloud and Dilated Cloud)
  var clearBit = 1 << 6;

  var isClear = qa.bitwiseAnd(clearBit).neq(0);   // bit 6 => eq 1
  var isFill = qa.bitwiseAnd(fillBit).neq(0);     // bit 0 => eq 1 => Fill
  var isCirrus = qa.bitwiseAnd(cirrusBit).neq(0); // bit 2 => eq 1 => Cirrus

  // Get pixels WHERE (Clear) AND (не Fill) AND (не Snow) AND (не Cirrus)
  var mask = isClear.and(isFill.not()).and(isCirrus.not());

  return image.updateMask(mask);
}

// Masking Landsat-5
var reflectanceColL5_nocloud = reflectanceColL5.map(maskCloudsL5);
// Masking Landsat-8
var reflectanceColL8_nocloud = reflectanceColL8.map(maskCloudsL8);

// Checking one image example
// var sampleMaskedL5 = reflectanceColL5_nocloud.filterMetadata('year', 'equals', example_yearL5).first();
// print('Example L5 no clouds', sampleMaskedL5);
// Map.addLayer(sampleMaskedL5, {}, 'L5 no clouds');

// var sampleMaskedL8 = reflectanceColL8_nocloud.filterMetadata('year', 'equals', example_yearL8).first();
// print('Example L8 no clouds', sampleMaskedL8);
// Map.addLayer(sampleMaskedL8, {}, 'L8 no clouds');

// Calculate NDVI
function getNDVI(image, nirband, redband) {
  var nir = image.select(nirband);
  var red = image.select(redband);

    var NDVI = image.expression(
    '(NIR - RED) / (NIR + RED)', {
      'NIR': nir,
      'RED': red
    }
  ).rename('NDVI');
    
    // Add NDVI band
    return image.addBands(NDVI);
}

// Apply function to reflectance VNIR+SWIR collection 
var NDVIcolL5 = reflectanceColL5_nocloud.map(function(image) {
  return getNDVI(image, "L_toa_B4", "L_toa_B3")
});
var NDVIcolL8 = reflectanceColL8_nocloud.map(function(image) {
  return getNDVI(image, "L_toa_B5", "L_toa_B4")
});

// Visualizing NDVI
// var sampleNDVI_L5 = NDVIcolL5.filterMetadata('year', 'equals', example_yearL5).first();
// Map.addLayer(sampleNDVI_L5, {bands: ['NDVI'], min: 0, max: 1}, 'NDVI L5');
//print('Processed collection NDVI L5:', NDVIcolL5);
// var sampleNDVI_L8 = NDVIcolL8.filterMetadata('year', 'equals', example_yearL8).first();
// Map.addLayer(sampleNDVI_L8, {bands: ['NDVI'], min: 0, max: 1}, 'NDVI L8');
//print('Processed collection NDVI L8:', NDVIcolL8);

///////////////////////////////////////////////////////////////////////////////////////

// Calculate NDBI and Built-up index BU

function addNDBI(image, nirband, swirband) {
  // Get channels from image
  var nir = image.select(nirband);
  var swir = image.select(swirband);
  // NDBI formula:
  // NDBI = ((SWIR - NIR) / (SWIR + NIR))
  var ndbi = image.expression(
    '(swir - nir) / (swir + nir)', {
      'swir': swir,
      'nir': nir
    }
  ).rename('NDBI');
  
  var ndvi = image.select("NDVI")
  // Calculate BU-index as ndbi-ndvi
  var bu = image.expression(
    'NDBI - NDVI', {
      'NDVI': ndvi,
      'NDBI': ndbi
    }
  ).rename('BU');
  
  // Add both indexes to image
  return image.addBands([ndbi, bu]);
}

// Apply function to NDVIcol collection
var NDVI_BU_colL5 = NDVIcolL5.map(function(image) {
  return addNDBI(image, "L_toa_B4", "L_toa_B5");
});
var NDVI_BU_colL8 = NDVIcolL8.map(function(image) {
  return addNDBI(image, "L_toa_B5", "L_toa_B6");
});

// Define the region for percentiles calculation
// and calculate FVC

var rgnperc = ee.FeatureCollection("projects/ee-olegkamalyuta/assets/mscpercentile").geometry();
Map.addLayer(rgnperc, {color: 'green'}, 'rgnpercentile');

////////////////////////////////////////////////////
// 1) Get NDVI maximum (multi-temporal composite)
////////////////////////////////////////////////////

// Selecting 'NDVI' band and get its absolute maximum over the time period
var maxNDVIImageL5 = NDVI_BU_colL5.select('NDVI').max();
var maxNDVIImageL8 = NDVI_BU_colL8.select('NDVI').max();

// 2) Calculating 10 and 90 percentiles for maxNDVIImage by region rgn
var percentiles_L5 = maxNDVIImageL5.reduceRegion({
  reducer: ee.Reducer.percentile([10, 90]),
  geometry: rgnperc,
  scale: 30,         // Landsat resolution
  bestEffort: true   // allows calculation simplifying for large regions
});
// same for Landsat 8
var percentiles_L8 = maxNDVIImageL8.reduceRegion({
  reducer: ee.Reducer.percentile([10, 90]),
  geometry: rgnperc,
  scale: 30,      
  bestEffort: true  
});

// Get calculated values
var NDVIsL5 = ee.Number(percentiles_L5.get('NDVI_p10'));
var NDVIvL5 = ee.Number(percentiles_L5.get('NDVI_p90'));
var NDVIsL8 = ee.Number(percentiles_L8.get('NDVI_p10'));
var NDVIvL8 = ee.Number(percentiles_L8.get('NDVI_p90'));

// Print them as numbers for verification
print('NDVIs (10th percentile) LS5:', NDVIsL5);
print('NDVIv (90th percentile) LS5:', NDVIvL5);
print('NDVIs (10th percentile) LS8:', NDVIsL8);
print('NDVIv (90th percentile) LS8:', NDVIvL8);

////////////////////////////////////////////////////
// 3) Calculate FVC for each image in NDVIcol
////////////////////////////////////////////////////
function addFVC(image, NDVIs, NDVIv) {
  // Get NDVI from image
  var ndvi = image.select('NDVI');
  // FVC formula:
  // FVC = ((NDVI - NDVIs) / (NDVIv - NDVIs))^2
  var fvc = ndvi.expression(
    '( (ndvi - ndviS) / (ndviV - ndviS) ) ** 2', {
      'ndvi': ndvi,
      'ndviS': NDVIs,
      'ndviV': NDVIv
    }
  ).rename('FVC');
  
  // If NDVI < NDVIs, set FVC = 0 --> classical approach
  fvc = fvc.where(ndvi.lt(NDVIs), 0);
  // If NDVI > NDVIv, set FVC = 1
  fvc = fvc.where(ndvi.gt(NDVIv), 1);
  
  // Add FVC to image
  return image.addBands(fvc);
}

// Apply function to NDVIcol collection
var FVCcolL5 = NDVI_BU_colL5.map(function(image) {
  return addFVC(image, NDVIsL5, NDVIvL5);
});
var FVCcolL8 = NDVI_BU_colL8.map(function(image) {
  return addFVC(image, NDVIsL8, NDVIvL8);
});

// Visualizing FVC
// var sampleImageFVC_L5 = FVCcolL5.filterMetadata('year', 'equals', example_yearL5).first();
// Map.addLayer(sampleImageFVC_L5, {bands: ['FVC'], min: 0, max: 1}, 'FVC L5');
// print('FVC collection:', FVCcolL5);
// var sampleImageFVC_L8 = FVCcolL8.filterMetadata('year', 'equals', example_yearL8).first();
// Map.addLayer(sampleImageFVC_L8, {bands: ['FVC'], min: 0, max: 1}, 'FVC L8');
// print('FVC collection:', FVCcolL8);

// Calculate epsilon (emissivity)
function addEmissivity(image) {
  // Set epsilon for bare soil and vegetation --> classical approach
  var epsilon_s = 0.97;
  var epsilon_v = 0.99;
  
  // Get FVC
  var fvc = image.select('FVC');
  
  // Calculate epsilon using formula eps_s * (1-fvc) + eps_v * fvc
  var epsilon = fvc.expression(
    'eps_s * (1 - fvc) + eps_v * fvc', {
      'eps_s': epsilon_s,
      'eps_v': epsilon_v,
      'fvc': fvc
    }
  ).rename('epsilon');
  
  // Adding epsilon band to an image
  return image.addBands(epsilon);
}

// Applying function to entire collection
var emissivityCol_L5 = FVCcolL5.map(addEmissivity);
var emissivityCol_L8 = FVCcolL8.map(addEmissivity);

// Visualizing emissivity
// var sampleImageEpsilon_L5 = emissivityCol_L5.filterMetadata('year', 'equals', example_yearL5).first();
// Map.addLayer(sampleImageEpsilon_L5, {bands: ['epsilon'], min: 0.97, max: 0.99}, 'Emissivity L5');
// print('Emissivity collection L5:', emissivityCol_L5);
// var sampleImageEpsilon_L8 = emissivityCol_L8.filterMetadata('year', 'equals', example_yearL8).first();
// Map.addLayer(sampleImageEpsilon_L8, {bands: ['epsilon'], min: 0.97, max: 0.99}, 'Emissivity L8');
// print('Emissivity collection L8:', emissivityCol_L8);

/////////////////////////////////////////////////////////////////////////////////////////
// Add water vapour content w to our images
/////////////////////////////////////////////////////////////////////////////////////////

// 1) Loading MERRA-2 reanalysis collection
var mERRA2 = ee.ImageCollection("NASA/GSFC/MERRA/slv/2")
  .select('TQV')
  .filterDate('1984-06-01', '2024-08-31');

// 2) Doing Join (saveFirst)
var timeFilter = ee.Filter.maxDifference({
  difference: 3600000,  // 1 hour in msec
  leftField: 'system:time_start',
  rightField: 'system:time_start'
});

var join = ee.Join.saveFirst({
  matchKey: 'merra'
});

var joinedFC_L5 = join.apply(emissivityCol_L5, mERRA2, timeFilter);
var joinedFC_L8 = join.apply(emissivityCol_L8, mERRA2, timeFilter);
// joinedFC is FeatureCollection, not ImageCollection!

// 3) Converting FeatureCollection to ImageCollection,
//    for each element get Landsat-image
//    and MERRA-image --> add w-band

function FeatureToImage(feature) {
  // feature — is ee.Feature, convert to ee.Image
    var landsatImg = ee.Image(feature);

    // Get MERRA-image
    var merraImg = ee.Image(feature.get('merra'));
    
    // Get TQV (vapour) and convert kg/m2 -> g/cm2
    var w = merraImg.select('TQV').multiply(0.1).rename('w');
    
    // Add w-band to Landsat-image
    return landsatImg.addBands(w);
}

var finalCollectionL5 = ee.ImageCollection(joinedFC_L5.map(FeatureToImage));
var finalCollectionL8 = ee.ImageCollection(joinedFC_L8.map(FeatureToImage));

// Now finalCollection is ImageCollection, where each image contains w-band with water vapour value

// Visualizing result
// var sampleL5 = finalCollectionL5.filterMetadata('year', 'equals', example_yearL5).first();
// print('Water Vapour image L5:', sampleL5);
// Map.addLayer(sampleL5, {bands: ['w'], min: 0, max: 3}, 'Water Vapour L5');
// var sampleL8 = finalCollectionL8.filterMetadata('year', 'equals', example_yearL8).first();
// print('Water Vapour image L8:', sampleL8);
// Map.addLayer(sampleL8, {bands: ['w'], min: 0, max: 3}, 'Water Vapour L8');


////////////// SAVE Cij coefficients for further ψ calculation /////////////////

// TIGR2311-based coefficients for Landsat-5 Band 6 by Jiménez-Muñoz
var scCoeffsL5 = {
  c11: 0.08158,   c12: -0.05707,  c13: 1.05991,
  c21: -0.58853,   c22: -1.08536,  c23: -0.00448,
  c31: -0.06201,  c32: 1.59086,   c33: -0.33513
};

// GAPRI4838-based coefficients for Landsat-8 Band 10 by Jiménez-Muñoz
var scCoeffsL8 = {
  c11: 0.04019,   c12:  0.02916,  c13: 1.01523,
  c21: -0.38333,   c22: -1.50294,  c23: 0.20324,
  c31: 0.00918,  c32: 1.36072,   c33: -0.27514
};

// Compute ψ1,2,3 parameters (psi1, psi2, psi3)
function addAtmosFunctions(image, scCoeffs) {
  // Get water vapour w
  var w = image.select('w');
  
  // Calculate psi1, psi2, psi3 using formula
  var psi1 = w.expression(
    'c11 * w * w + c12 * w + c13', {
      'w': w,
      'c11': scCoeffs.c11, 'c12': scCoeffs.c12, 'c13': scCoeffs.c13
    }
  ).rename('psi1');
  
  var psi2 = w.expression(
    'c21 * w * w + c22 * w + c23', {
      'w': w,
      'c21': scCoeffs.c21, 'c22': scCoeffs.c22, 'c23': scCoeffs.c23
    }
  ).rename('psi2');
  
  var psi3 = w.expression(
    'c31 * w * w + c32 * w + c33', {
      'w': w,
      'c31': scCoeffs.c31, 'c32': scCoeffs.c32, 'c33': scCoeffs.c33
    }
  ).rename('psi3');
  
  // Add bands psi1, psi2, psi3 to image
  return image.addBands([psi1, psi2, psi3]);
}

// Apply function addAtmosFunctions for every image in collection
var withPsiL5 = finalCollectionL5.map(function(image) {
  return addAtmosFunctions(image, scCoeffsL5);
});

var withPsiL8 = finalCollectionL8.map(function(image) {
  return addAtmosFunctions(image, scCoeffsL8);
});

// Visualizing ψ1,2,3
//var sampleImgpsiL5 = withPsiL5.filterMetadata('year', 'equals', example_yearL5).first();
//Map.addLayer(sampleImgpsiL5, {bands: ['psi1', 'psi2', 'psi3'], min: 0, max: 2}, 'psi L5');

//var sampleImgpsiL8 = withPsiL8.filterMetadata('year', 'equals', example_yearL8).first();
//Map.addLayer(sampleImgpsiL8, {bands: ['psi1', 'psi2', 'psi3'], min: 0, max: 2}, 'psi L8');

// Print results for verifying
// print('Psi1,2,3 L5:', sampleImgpsiL5);
// print('Psi1,2,3 L8:', sampleImgpsiL8);

//// Prepare γ and δ variables for LST calculation:

var byL5 = 1256 // = c2/λ_effective = const for sensor
var byL8 = 1324

function getvarsdY(image, by) {
  
  // Get Tb and L lambda from image for further calculations
  var Tsen = image.select('Tb_TIR');
  var Lsen = image.select('L_lambda_TIR')
  
  // Calculate d and Y using Jiménez-Muñoz formulas
  var d = image.expression(
    'Tsen - (Tsen * Tsen / by)', {
      'Tsen': Tsen,
      'by': by
    }
  ).rename('d');
  
  var Y = image.expression(
    '(Tsen * Tsen) / (by * Lsen)', {
      'Tsen': Tsen,
      'by': by,
      'Lsen': Lsen
    }
  ).rename('Y');

  
  // Add variables to image
  return image.addBands([d,Y]);
}

// Add variables to image containing ψ1,2,3
var withdY_L5 = withPsiL5.map(function(image) {
  return getvarsdY(image, byL5);
});

var withdY_L8 = withPsiL8.map(function(image) {
  return getvarsdY(image, byL8);
});

// Visualising d and Y
// var sampleImgdY_L5 = withdY_L5.filterMetadata('year', 'equals', example_yearL5).first();
// Map.addLayer(sampleImgdY_L5, {bands: ['d', 'Y']}, 'd, Y L5');
// print('vars d, Y L5', sampleImgdY_L5);
// var sampleImgdY_L8 = withdY_L8.filterMetadata('year', 'equals', example_yearL8).first();
// Map.addLayer(sampleImgdY_L8, {bands: ['d', 'Y']}, 'd, Y L8');
// print('vars d, Y L8', sampleImgdY_L8);

//////////////////////////////////////////////////////////////////////////////////
///////////FINAL LST CALCULATION using variables computed earlier////////////////
//////////////////////////////////////////////////////////////////////////////////

function getLST(image) {
  
  var Y = image.select('Y');
  var epsilon = image.select('epsilon');
  var psi1 = image.select('psi1');
  var psi2 = image.select('psi2');
  var psi3 = image.select('psi3');
  var Ltoa = image.select('L_lambda_TIR');
  var d = image.select('d');
  
  var LST = image.expression(
      'Y * ((1 / epsilon) * (psi1 * Ltoa + psi2) + psi3) + d',
      {
        'Y': Y,
        'epsilon': epsilon,
        'psi1': psi1,
        'psi2': psi2,
        'psi3': psi3,
        'Ltoa': Ltoa,
        'd': d
      }
      ).rename('LST');
  
// Add LST to image
  return image.addBands([LST]);
}

// Add LST to image
var LST_L5 = withdY_L5.map(getLST);
var LST_L8 = withdY_L8.map(getLST);

// Visualizing LST 
// var sampleImgLST_L5 = LST_L5.filterMetadata('year', 'equals', example_yearL5).first();
// Map.addLayer(sampleImgLST_L5, {bands: ['LST']}, 'LST L5 example');
// print('LST L5 example', sampleImgLST_L5);
// var sampleImgLST_L8 = LST_L8.filterMetadata('year', 'equals', example_yearL8).first();
// Map.addLayer(sampleImgLST_L8, {bands: ['LST']}, 'LST L8 example');
// print('LST L8 example', sampleImgLST_L8);

////////////////////////////Merge collection///////////////////////////////////////////

var merged_example_year = 2022 // Visualization year for merged collection of LS5,8

// Merge L5 and L8 collections
var mergedLSTCollection = ee.ImageCollection(LST_L5.merge(LST_L8));

// Verifying size of general collection
// print('Merged LST collection (L5 + L8):', mergedLSTCollection);
print('Total images number in Merged Collection:', mergedLSTCollection.size());

var sampleImgLSTMerged = mergedLSTCollection.filterMetadata('year', 'equals', merged_example_year).first();
Map.addLayer(sampleImgLSTMerged, {bands: ['LST'], min: 273, max: 320}, 'General merged LST example');
// print('Merged collection example:', sampleImgLSTMerged);

/////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Aggregation beginning /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

var roi = ee.FeatureCollection("projects/ee-olegkamalyuta/assets/ckad").geometry();
Map.addLayer(roi, {color: 'blue'}, 'roi_rgn');

// Group collection by year and then calculate mean LST by year

// Aggregate images by year
function aggregateByYear(collection, year) {
  // Filter collection by year
  var yearlyCollection = collection.filterMetadata('year', 'equals', year);
  
  // Check if images exist in chosen year
  var hasImages = yearlyCollection.size().gt(0);
  
  // Returning mean image only if data exist using ee.Algorithms.If
  return ee.Algorithms.If(
    hasImages,
    yearlyCollection.map(function(image) {
      return image.select(['LST', 'NDBI', 'NDVI', 'BU']).clip(roi)
    }).mean().set('year', year),
    null // Return null if there is no data for chosen year
  );
}

// Unique years list
var years = ee.List.sequence(L5_START_YEAR, L8_END_YEAR);

// Apply aggregation for each year
var aggregatedLSTCollection = years.map(function(year) {
  return aggregateByYear(mergedLSTCollection, year);
});

// Convert list to ImageCollection
aggregatedLSTCollection = ee.ImageCollection(aggregatedLSTCollection);

// Verify number of images in aggregated collection
print('Number of images aggregated by year:', aggregatedLSTCollection.size());

// Visualizing results
var sampleImage = aggregatedLSTCollection.filterMetadata('year', 'equals', merged_example_year).first();
Map.addLayer(sampleImage, {bands: ['LST'], min: 273, max: 320}, 'Aggregated LST example');
print('Aggregated LST Image example', sampleImage);


// !!! LST images have been received !!!

////////// Analysis part: compute mean LST by roi region based on aggregated for every year data 
//////////                for further validation with reanalysis

function getMeanForRegion(image) {
  var lst = image.select('LST');
  var year = image.get('year');
  
  var meanRegion = image.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: roi, // ROI
      scale: 30,
      maxPixels: 1e10
    });
    
  var meanLSTbyRegion = ee.Algorithms.If(meanRegion.contains('LST'), meanRegion.get('LST'), null);
  
  return ee.Feature(null, {
      'year': year,
      'mean_LST_by_Region': meanLSTbyRegion
    });
  
}

var meanLSTbyRegion = aggregatedLSTCollection.map(getMeanForRegion);


// Export CSV table with ROI-averaged LST value for every year
Export.table.toDrive({
  collection: meanLSTbyRegion,
  description: 'LST_Mean_by_ROI',
  fileFormat: 'CSV',
  folder: 'kursach_lst_poly', // your path
});

print("Average LST across the entire ROI:", meanLSTbyRegion.limit(10));

/////////////////////////////////// Calculate number of aggregated LST /////////////////////////////////////////

// Count images number for general aggregated collection
// var yearCountsLST = countImagesByYear(aggregatedLSTCollection, years);

// print("Number of aggregated by year LST:", yearCountsLST);
// Ideally should be 1 image per year but may be some gaps

var yearsWithData = aggregatedLSTCollection.aggregate_array('year').sort();
var yearsWithoutData = years.removeAll(yearsWithData);
print('Gap years in Final Aggregated LST Collection:', yearsWithoutData);

/////////////
///////////// Analysis: normalize EVERY LST-image using mean rural LST to get SUHI values
/////////////


// Calculate SUHI
function normalizeAllImage(image, rural) {
  // Get LST band from image clipped by ROI 
  var lstBand = image.select('LST').clip(roi);
  
  // Compute mean LST in rural area
  var ruralMean = lstBand.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: rural,
    scale: 30,
    maxPixels: 1e10
  }).get('LST');
  
  // LST normalizing: subtraction mean rural LST for a year
  var normalizedLST = lstBand.subtract(ee.Image.constant(ruralMean)).toFloat();
  
  // Get BU and NDVI-band (original)
  var buBand = image.select('BU').toFloat();
  var ndviBand = image.select('NDVI').toFloat();
  
  // Add BU-band to image
  var normalizedImage = normalizedLST.addBands([buBand, ndviBand]).set('year', image.get('year'));
  
  return normalizedImage;
}

var SUHICollection = aggregatedLSTCollection.map(function(image)
{
  return normalizeAllImage(image, rural)
}
)

print('Size of COMPLETE collection of normalized (SUHI) images:', SUHICollection.size());

var sample_normalized_one = SUHICollection.first();
Map.addLayer(sample_normalized_one, {bands: ['LST']}, 'SUHI example');

// print('COMPLETE collection of normalized (SUHI) images:', SUHICollection)


//// Calculate SUHI areas for each agglomeration part 
//// (T>2 threshold used taking into account the error of the LST obtaining algorithm)

// Load polygons of parts of agglomeration
var agglomeration = ee.FeatureCollection("projects/ee-olegkamalyuta/assets/agglomeration");
Map.addLayer(agglomeration, {color: 'blue'}, 'agglomeration');

// Calculate area of LST > 2 for image and agglomeration polygon
function calculateArea(image, polygon) {
  var year = image.get('year');
  var agl = polygon.get('agl');
  
  // Clip image by the polygon of the agglomeration part
  var clippedImage = image.clip(polygon);
  
  // Get pixels where LST > 2
  var hotPixels = clippedImage.select('LST').gt(2);
  
  // Get area in km² in appropriate projection
  var area = hotPixels.multiply(ee.Image.pixelArea()).reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: polygon.geometry(),
    scale: 30,
    crs: 'EPSG:32637',  // UTM 37N is suitable for Moscow
    maxPixels: 1e10
  }).get('LST');
  
  // Return feature containing area value
  return ee.Feature(null, {
    year: year,
    agl: agl,
    square: ee.Number(area).divide(1e6)  // Convert m² to km²
  });
}

// Get polygons list
var polygons = agglomeration.toList(agglomeration.size());

// Calculate area for every image
function processImage(image) {
  var areas = polygons.map(function(polygon) {
    return calculateArea(image, ee.Feature(polygon));
  });
  return ee.FeatureCollection(areas);
}

// Flattening collection to get appropriate results
var results = SUHICollection.map(processImage).flatten();

// Export area by part of agglomeration in CSV
Export.table.toDrive({
  collection: results,
  description: 'SUHI_Areas',
  folder: 'kursach_lst_poly', // your path
  fileFormat: 'CSV',
  selectors: ['year', 'agl', 'square']
});

//////////////////////////////////////////////////////////////////////////////////////////
// Analysis: beginning of processing by Local Climate Zones LCZ polygons
//////////////////////////////////////////////////////////////////////////////////////////

// Load LCZs
var lcz = ee.FeatureCollection("projects/ee-olegkamalyuta/assets/lcz_zones");
Map.addLayer(lcz, {color: 'pink'}, 'lczs');

// Compute mean LST and indexes for LCZ polygon and year
function computeMeanLST(feature) {
  var polygon = feature.geometry();
  var lcz_id = feature.get('id');
  var lcz = feature.get('lcz')
  var lcz_name = feature.get('name');
  var rus_name = feature.get('rus')

  // Get mean for each image (year) for current polygon
  var yearlyMeans = aggregatedLSTCollection.map(function(image) {
    var year = image.get('year');

    // Clip image by polygon and compute mean
    var meanDict = image.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: polygon,
      scale: 30,
      maxPixels: 1e8
    });

    // Get mean band values using dict keys
    var meanLST = meanDict.get('LST');
    var meanNDBI = meanDict.get('NDBI');
    var meanBU = meanDict.get('BU');
    
//    var meanLST = ee.Algorithms.If(meanDict.contains('LST'), meanDict.get('LST'), null);
//    var meanNDBI = ee.Algorithms.If(meanDict.contains('NDBI'), meanDict.get('NDBI'), null);
//    var meanBU = ee.Algorithms.If(meanDict.contains('BU'), meanDict.get('BU'), null);

    return ee.Feature(null, {
      'year': year,
      'lcz_id': lcz_id,
      'lcz': lcz,
      'lcz_name': lcz_name,
      'rus_name': rus_name,
      'mean_LST': meanLST,
      'mean_NDBI': meanNDBI,
      'mean_BU': meanBU
    });
  });

  return yearlyMeans;
}

// Apply for each polygon and flatten results
var finalResult = lcz.map(computeMeanLST).flatten();

// Export mean LST for each LCZ by year in CSV
Export.table.toDrive({
  collection: finalResult,
  description: 'LST_Mean_by_LCZ',
  fileFormat: 'CSV',
  folder: 'kursach_lst_poly', // your path
});

// Print first lines for verifying
print('Average LST by LCZ examples:', finalResult.limit(10));

/////////////// Visualization: Create and export images collection averaged over 5 years ////////////////

// Define 5-year periods
var periods = [
  {start: 1984, end: 1988},
  {start: 1989, end: 1993},
  {start: 1994, end: 1998},
  {start: 1999, end: 2003},
  {start: 2004, end: 2008},
  {start: 2009, end: 2014},
  {start: 2015, end: 2019},
  {start: 2020, end: 2024}
];

// Calculate averaged image
function createMeanImage(period) {
  var startYear = ee.Number(period.start);
  var endYear = ee.Number(period.end);
  
  // Filter collection by periods
  var filteredCollection = aggregatedLSTCollection.filter(
    ee.Filter.rangeContains('year', startYear, endYear)
  );
  
  // Compute mean LST, BU and NDVI
  var meanImage = filteredCollection.select('LST').mean();
  var meanImageBU = filteredCollection.select('BU').mean();
  var meanImageNDVI = filteredCollection.select('NDVI').mean();
  
  // Get string for period property in the format "startYear_endYear"
  var periodString = startYear.format('%d').cat('_').cat(endYear.format('%d'));
  
  // Add period property for period identifying
  meanImage = meanImage.addBands([meanImageBU, meanImageNDVI]).set('period', periodString);
  
  return meanImage;
}

// Create list of averaged images for each period
var meanImages = periods.map(createMeanImage);

// Convert list to ImageCollection
var meanImagesCollection = ee.ImageCollection(meanImages);

// Export every image to Google Drive
meanImagesCollection.getInfo(function(info) {
  var images = info.features;
  images.forEach(function(imageInfo) {
    var image = ee.Image(meanImagesCollection.filter(
      ee.Filter.eq('period', imageInfo.properties.period)
    ).first());
    var period = imageInfo.properties.period;
  
    // Настройка экспорта
    Export.image.toDrive({
      image: image,
      description: 'LST_BU_NDVI_Mean_' + period,      
      folder: 'kursach_lst_poly',   // your path           
      fileNamePrefix: 'LST_BU_NDVI_Mean_' + period,    // Prefix
      region: roi,                            // ROI
      scale: 30,                               // Landsat resolution
      crs: 'EPSG:4326',                        // CRS WGS84
      maxPixels: 1e10                          
    });
  });
});

// Print collection for verifying 
print('Aggregated images collection:', meanImagesCollection);

//////////////////////// Analysis: Calculate SUHI aggregated over 5-year periods //////////////////

var rural = ee.FeatureCollection("projects/ee-olegkamalyuta/assets/rural").geometry();
Map.addLayer(rural, {color: 'green'}, 'rural_rgn');

// New normalization function
function normalizeImage(image, rural) {
  // Get LST band
  var lstBand = image.select('LST');
  
  // Compute mean LST in Rural area
  var ruralMean = lstBand.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: rural,
    scale: 30,
    maxPixels: 1e10
  }).get('LST');
  
  // Normalize LST again
  var normalizedLST = lstBand.subtract(ee.Image.constant(ruralMean)).toFloat();
  
  // Get BU-band
  var buBand = image.select('BU').toFloat();
  var ndviBand = image.select('NDVI').toFloat();
  
  // Get SUHI image for every period
  var normalizedImage = normalizedLST.addBands([buBand, ndviBand]).set('period', image.get('period'));
  
  return normalizedImage;
}

// Apply normalization
var normalizedImagesCollection = meanImagesCollection.map(function(image) {
  return normalizeImage(image, rural);
});

// Export normalized images to Google Drive
normalizedImagesCollection.getInfo(function(info) {
  var images = info.features;
  images.forEach(function(imageInfo) {
    var image = ee.Image(normalizedImagesCollection.filter(
      ee.Filter.eq('period', imageInfo.properties.period)
    ).first());
    var period = imageInfo.properties.period;
    
    Export.image.toDrive({
      image: image,
      description: 'LST_BU_NDVI_Mean_Normalized_' + period,
      folder: 'kursach_lst_poly',   // your path
      fileNamePrefix: 'LST_BU_NDVI_Mean_Normalized_' + period,
      region: roi, // ROI
      scale: 30,    // Landsat resolution
      crs: 'EPSG:4326', // CRS WGS84
      maxPixels: 1e10
    });
  });
});

// Visualization for normalized images
var sampleImage = normalizedImagesCollection.first();
Map.addLayer(sampleImage, {bands: ['LST'], min: -10, max: 25}, 'Normalized Aggregated LST');
print('Normalized Aggregated LST example', sampleImage);

// Print collection for verifying
print('Normalized Aggregated LST collection:', normalizedImagesCollection);


///////////////////////////////////// Visualization and analysis: clustering /////////////
//////////////////////////////////////////////////////////////////////////////////////////

// Training Data generation

// function clusterImage(image) {
  
//   // var trainingData = image.sample({ // Random data without specifying initial points
//   //   region: roi,    // ROI
//   //   scale: 30,       // Landsat resolution
//   //   numPixels: 10e4, // Pixels number
//   //   seed: 0          // Seed
//   // });
  
//   var trainingData = image.sample({ // Clipping using lcz --> clusterizer is trained on lcz zones pixels
//     region: lcz,
//     scale: 30,       // Landsat resolution
//     numPixels: 10e4, // Pixels number
//     seed: 0          // Seed
//   });
  
//   var clusterer = ee.Clusterer.wekaKMeans(4).train({
//   features: trainingData,
//   inputProperties: ['LST']  // Clustering band
//   });
  
//   var clusteredImage = image.select('LST').cluster(clusterer).rename('clusters')
  
  
//   return clusteredImage.set('period', image.get('period'));
// }

// // Clustering
// var clusteredImagesCollection = normalizedImagesCollection.map(clusterImage)

// // Visualizing clustering results example
// var sampleClusterised = clusteredImagesCollection.first()

// Map.addLayer(sampleClusterised, {
//     bands: ['clusters'],
//     min: 0,
//     max: 3,
//     palette: ['#a9ffac', '#ff0000', '#4fb7d0', '#ffcf45']
//   }, 'SUHI clusters example');
  
// print('Clusterized images number:', clusteredImagesCollection.size())
// // print('Clusterized aggregated images collection:', clusteredImagesCollection)

// // Export clustering results
// clusteredImagesCollection.getInfo(function(info) {
//   var images = info.features;
//   images.forEach(function(imageInfo) {
//     var image = ee.Image(clusteredImagesCollection.filter(
//       ee.Filter.eq('period', imageInfo.properties.period)
//     ).first());
//     var period = imageInfo.properties.period;
    
//     Export.image.toDrive({
//       image: image,
//       description: 'Clusters_lcz_Normalized_' + period,
//       folder: 'kursach_lst_poly',   // your path
//       fileNamePrefix: 'Clusters_lcz_Normalized_' + period,
//       region: roi,   // ROI
//       scale: 30,      // Landsat resolution
//       crs: 'EPSG:4326', // CRS WGS84
//       maxPixels: 1e10
//     });
//   });
// });
