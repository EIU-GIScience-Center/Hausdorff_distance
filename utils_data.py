# -*- coding: utf-8 -*-
"""
Created on Sat May 25 19:21:05 2024

@author: Administrator
"""

import shapefile
import win32clipboard

def get_sample_feature(feat_name):
    feat_name = feat_name.lower()
    def feature_from_data_folder(sf_name):
        # folder = r"C:\CaGIS Board Dropbox\cantaloupe bob\Barry\Research\Projects\polyline difference metrics\Hausdorff\data\Shape Files"
        folder = r"sample_data\original"
        filename = f"{folder}\\{sf_name}.shp"
        sf = shapefile.Reader(filename)
        return sf.shapes()[0].points
    if feat_name == "cannonball":
        return feature_from_data_folder("Cannonball")
    elif feat_name == "india_bangladesh":
        return feature_from_data_folder("India_Bangladesh")
    elif feat_name.lower() == "lake_shelbyville":
        return feature_from_data_folder("lake_shelbyville")
    elif feat_name == "manhattan":
        return feature_from_data_folder("Manhattan")

def copy2clip(txt):
    # set clipboard data
    win32clipboard.OpenClipboard()
    win32clipboard.EmptyClipboard()
    win32clipboard.SetClipboardText(txt)
    win32clipboard.CloseClipboard()

def create_polyline_shapefile(polylines,filepath):
    """
    Writes geometry data representing a set of possibly multipart polylines 
    into a shapefile 

    Parameters
    ----------
    polylines : list of lists of lists of (x,y) tuples
        The geometry of the polylines.
    filepath : str
        The full path of the file to save to.

    Returns
    -------
    None.

    """
    # Create a point shapefile
    w = shapefile.Writer(filepath, shapeType=shapefile.POLYLINE)
    
    # Define fields
    w.field('ID', 'N') 
    
    for i in range(len(polylines)):
        # Add polylines and records""
        w.line(polylines[i])
        w.record(i)

    # Save the shapefile
    w.close()