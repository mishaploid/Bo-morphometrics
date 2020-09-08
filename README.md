# Brassica oleracea leaf morphometrics

Analysis of leaf shape for a diversity panel of Brassica oleracea and wild C genome relatives to examine differences across developmental stages (leaves 2, 3, and 4) and accessions. 

## Getting Started

Recommend following the excellent guidleines provided by Laura Klein and Harlan Svoboda (https://bio-protocol.org/e2269)

### Prerequisites

All open source software 
ImageJ (https://imagej.nih.gov/ij/)
SHAPE (http://lbm.ab.a.u-tokyo.ac.jp/iwata/shape/index.html) - requires Winebottler if not installing on Windows (http://winebottler.kronenberg.org/)
R (https://www.r-project.org/)

```
Examples 

```
### Convert scans to binary images and crop individual leaves

ImageJ macro 

## Convert TIF binary images to BMP with RGB

Image J > Process > Macro

Select Input/Output directories 

Paste text from src/convert_to_bmp.txt


### Rename files to specify leaf number based on position from ImageJ

src/edit_filenames.Rmd (R Notebook) 

src/rename_files.sh (bash script to rename files based on output from edit_filenames.Rmd)

### Landmarking 

![Landmark positions](reports/Landmarking.png)


## Authors

* **Sarah Turner-Hissong**
* **Makenzie Mabry**
* **Evan Gallagher**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* https://bio-protocol.org/e2269
* Dan Chitwood

 
