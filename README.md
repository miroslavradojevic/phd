# Automated neuron tracing using probability hypothesis density filtering #
##### ImageJ plugin #####

Radojević, M., & Meijering, E. (2017). Automated neuron tracing using probability hypothesis density filtering. Bioinformatics, 33(7), 1073-1080.

[doi.org/10.1093/bioinformatics/btw751](https://doi.org/10.1093/bioinformatics/btw751)

### Usage ###
* copy [phd_.jar](https://github.com/miroslavradojevic/phd/releases/download/1.0.0/phd_1.0.0.jar) and [vess_.jar](https://github.com/miroslavradojevic/vess/releases/download/1.0.0/vess_1.0.0.jar) to ImageJ's plugins directory
* Open ImageJ > Plugins > BrainCadet > PHD
* select file (works with 8bit images and image stacks, dark background)
* enter values in the parameter menu

### ImageJ macro for batch processing ###
* parameter grid 
* possibility to use the parameter grid to process all existing *.tif 8bit images in given directory 
* ij macro: ij-macro/run_phd.ijm
* test case (~/test/ directory with run_phd.ijm and img.tif), example call: 
```
#!java

java -Xmx4g -jar ~/ImageJ/ij.jar -ijpath ~/ImageJ/plugins/ -batch ~/test/run_phd.ijm ~/test/
```
