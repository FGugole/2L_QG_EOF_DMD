# 2 layers QG model with stochastic forcing 
This project incudes the implementation of a numerical model solving the 2-layer Quasi-Geostrophic (QG) model 
with stochastic forcing as described in this [paper](https://www.essoar.org/doi/abs/10.1002/essoar.10501115.1). 
For details on the model and its numerical discretization/implementation please see the aforementioned paper. 

## Usage
This repository includes three main files: QG_det.f90, QG_sto.f90 and QG_sto_DMDpropagation.f90. 
All of them make use of the **FFTW** module (which is not provided here but can be easily found at www.fftw.org), 
hence be sure to link it correctly. QG_sto.f90 and QG_sto_DMDpropagation.f90 require also the LAPACK library.

**QG_det.f90** corresponds to the discretization of the **deterministic** 2-layer QG model, hence it does not 
include any stochastic forcing. 
>QG_det.f90 requires to be linked to the modules/subroutines: AJ4.f90, EnergyOutput.f90, IO_restart.f90 and Spectra.f90.

**QG_sto.f90** includes the **stochastic** forcing but does not allow for the propagation of the DMD modes. 
**QG_sto_DMDpropagation.f90** includes the stochastic forcing and allows for the propagation of the DMD modes. 
Both QG_sto.f90 and QG_sto_DMDpropagation.f90 recovers the results of QG_det.f90 when the stochastic forcing 
is set to zero, but have higher computational times. 
>Both QG_sto.f90 and QG_sto_DMDpropagation.f90 require to be linked to the modules/subroutines: AJ4.f90, 
BoxMuller.f90, DMD.f90, EnergyOutput.f90, IO_restart.f90, mtfort.f90, NoiseStructure.f90, Projection.f90 and Spectra.f90.

## Project status
The development of the project has been currently slowed down. If interested in maintainig actively this project, 
feel free to get in contact with me: federica.gugole@uni-hamburg.de

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)
