import rasterio
from rasterio.io import MemoryFile
from rasterio.mask import mask
import numpy as np
from shapely.geometry import box
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.crs import CRS


              
# Load Sentinel-2 bands (Red: B4, NIR: B8) from JP2 files
def load_bands(red_path, nir_path):
    with rasterio.open(red_path) as red_src:
        red = red_src.read(1)
        red_meta = red_src.meta
    
    with rasterio.open(nir_path) as nir_src:
        nir = nir_src.read(1)
        nir_meta = nir_src.meta
    
    return red-1000, nir-1000, red_meta

# Calculate NDVI
def calculate_ndvi(red, nir):
    ndvi = (nir.astype(float) - red.astype(float)) / (nir + red)
    ndvi = np.clip(ndvi, -1, 1)  # Clip values to range [-1, 1]
    return ndvi

# Clip the NDVI raster to a smaller extent (bounding box)
def clip_ndvi(ndvi, ndvi_meta, bbox):
    # Create a MemoryFile to store the NDVI raster in memory
    with MemoryFile() as memfile:
        # Write the NDVI array into the MemoryFile as a raster
        with memfile.open(**ndvi_meta) as ndvi_src:
            ndvi_src.write(ndvi, 1)
        
        # Now, reopen the MemoryFile in read mode and apply the mask
        with memfile.open() as ndvi_src:
            geometry = [box(*bbox)]
            clipped_ndvi, clipped_transform = mask(ndvi_src, geometry, crop=True)
        
        # Update metadata to reflect the new size and transform
        clipped_meta = ndvi_meta.copy()
        clipped_meta.update({
            "height": clipped_ndvi.shape[1],
            "width": clipped_ndvi.shape[2],
            "transform": clipped_transform
        })
    
    return clipped_ndvi, clipped_meta

# Save raster as a GeoTIFF
def save_geotiff(output_path, data, meta):
    meta.update(driver='GTiff')  # Set driver to GeoTIFF
    with rasterio.open(output_path, 'w', **meta) as dst:
        dst.write(data, 1)

# Main workflow
def process_sentinel_scene(red_path, nir_path, output_path, bbox):
    # Load the red and NIR bands
    red, nir, red_meta = load_bands(red_path, nir_path)
    
    # Calculate NDVI
    ndvi = calculate_ndvi(red, nir)
    
    # Create metadata for NDVI raster (single band)
    ndvi_meta = red_meta.copy()
    ndvi_meta.update({"count": 1, "dtype": 'float32', "nodata": None, "driver": "GTiff"})  # NDVI is single-band, and in float32 format
    # ndvi_meta.update(driver='GTiff')
    # Clip NDVI to the specified extent
    clipped_ndvi, clipped_meta = clip_ndvi(ndvi, ndvi_meta, bbox)
    
    # Save the clipped NDVI as a new GeoTIFF
    save_geotiff(output_path, clipped_ndvi[0,:,:], clipped_meta)
    
    with rasterio.open(output_path, 'r+') as dataset:
        # Assign the new CRS to the raster
        dataset.crs = CRS.from_epsg(25832)  

    
# Define file paths, bounding box, and target EPSG
f1='S2B_MSIL2A_20240720T103629_N0510_R008_T32VNH_20240720T121155.SAFE'
f2='L2A_T32VNH_A038501_20240720T103633'
f3_4='T32VNH_20240720T103629_B04_10m'
f3_8='T32VNH_20240720T103629_B08_10m'

red_band_path = f1+'/GRANULE/'+f2+'/IMG_DATA/R10m/'+f3_4+'.jp2'  # Sentinel-2 Red band (B4) in JP2 format
nir_band_path = f1+'/GRANULE/'+f2+'/IMG_DATA/R10m/'+f3_8+'.jp2'  # Sentinel-2 NIR band (B8) in JP2 format
output_path = f1+'/output_ndvi_20240720.tif'


bbox = [549100, 6267800, 574200, 6287800]  # Define the bounding box for clipping

# Process the Sentinel-2 scene
process_sentinel_scene(red_band_path, nir_band_path, output_path, bbox)
