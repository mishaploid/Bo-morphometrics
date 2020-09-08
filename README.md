# Morphometric analysis of Brassica oleracea leaves 

This repo includes scripts and processed data for a morphometric analysis of leaf scans in a diversity panel of Brassica oleracea and wild C genome relatives. The goal is to examine variation in leaf shape across developmental stages (leaves 2, 3, and 4) and accessions. 

## Prerequisites
* R (https://www.r-project.org/)

```
Examples 

```
## Image preprocessing steps 

1. Convert scans to binary images and crop individual leaves (python/PlantCV)
2. Export traditional morphometric descriptors (R/Momocs)
3. Compute Elliptical Fourier Descriptors (EFDs) to describe leaf contours (R/Momocs)
4. Landmark images to facilitate alignment 
   <img src="reports/Landmarking.png" width="200" height="200"> 


## Authors

* **Sarah Turner-Hissong**
* **Makenzie Mabry**
* **Evan Gallagher**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Resources
[Comprehensive Methods for Leaf Geometric Morphometric Analyses](https://bio-protocol.org/e2269) by Laura Klein and Harlan Svoboda 


## Acknowledgments

* https://bio-protocol.org/e2269
* Dan Chitwood

 
