# Polyline Hausdorff 
Pure python code to compute the Hausdorff distance between two polylines.

# What is the Hausdorff distance between two polylines?
The Hausdorff distance between polylines is the length of the largest gap between them. It may be expressed more formally as the maximum minimum distance, i.e. the maximum of **d(a&#8594;B)** and **d(b&#8594;A)** for all points **a** and **b** on polylines **A** and **B** respectively. 

# How does this code differ from SciPy/Shapely?
The SciPy and Shapely libraries do not compute the true Hausdorff distances between polylines.
* SciPy's [directed Hausdorff distance function](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.directed_hausdorff.html) computes the Hausdorff distance between two point sets. If given two polylines, it will compute the maximum distance between a **vertex** on one polyline and the nearest **vertex** on the other polyline. 
* Shapely's [Hausdorff distance](https://shapely.readthedocs.io/en/stable/manual.html#general-attributes-and-methods) function is not well documented, but [our experiments](https://generalisation.icaci.org/downloads/abs2019/Abs2019_paper_8.pdf) have shown that it actually calculates the maximum minimum distance between a ***vertex*** on one polyline and ***any point*** on the other polyline.
* In contrast, the Polyline Hausdorff function computes the true Hausdorff distance, which is the maximum minimum distance between ***any point*** on one polyline and ***any point*** on the other polyline.

The following figure illustrates the difference between the SciPy, Shapely and Polyline Hausdorff distances:

![Hausdorff comparison](https://user-images.githubusercontent.com/13248690/121942432-15951b00-cd16-11eb-9d6c-ea787f154f66.png)

# How do I use the code?
Basic usage is illustrated in the file, `illustration.py`.

# How does it work? 
The code follows the general algorithm described in:

[J.F. Hangou&#0235;t (1995). Computation of the Hausdorff distance between plane vector polylines. *Proceedings of the International Symposium on Computer-Assistated Cartography (AutoCarto XII).*](https://cartogis.org/docs/proceedings/archive/auto-carto-12/pdf/computation-of-the-hausdorff-distance-between-plane.pdf)

The algorithm has been updated in several ways to ensure robustness and improve computational efficiency. Details on these improvements are currently being prepared for publication.

# Current Status (updated Aug: 2024)
This is still a work in process. The code is fully functional and has been tested against a variety of test cases. The code still requires further cleaning to remove legacy functions bot main code is in a mature state. This code comes without any warranty, and is recommended for research purposes only.

## Robustness
The code has been tested on several cases including coincident, overlapping and perpendicular line segments with random transformations to catch problems involving floating point precision arithmetic. It appears to be robust but further testing is warranted for production-level usage. 

## Computational Efficiency
Computational efficiency has been improved with a "three-circle-filter" and is usually subquadratic but approaches theoretical cubic running time in the very unlikely worst case scenario. 

## Code Testing
Test functions need to be modified to reflect refactoring of main code.

# Who created this?
This code is being developed by Dr. Barry Kronenfeld at Eastern Illinois University with support from the USGS. However this repository has not been approved or endorsed by EIU or the USGS. 

Much of the code for the Hausdorff distance calculation was developed by students in the Spring 2021 section of **GEO 4910: GIS Programming** at Eastern Illinois University. Team members are:

* Luke Jansen
* Tanner Jones
* Farouk Olaitan
* Megshi Thakurtimes

# Is this free to use?
Yes, it is licensed with the open source MIT License. If you find this useful, acknowledgement would be appreciated.
