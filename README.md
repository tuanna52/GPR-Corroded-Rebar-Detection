# Smart System for Non-Invasive Inspection of Corroded Metal Structures in Bridges
This project aims at creating an electromagnetic model from first principles aiming at detecting corroded and pristine metal structures embedded into concrete and asphalt. The simulated structures of bridges in the experiment were examined by Frequency Modulated Continuous Wave (FMCW) radar. The received signal was processed by using model-based technique, which based on both the prior knowledge of how electromagnetic waves behave when reaching different interfaces and curve-fitting algorithms. This model-based approach showed promising results in studying the underlying physical properties of each layer below the surface of the bridges.

# Files' Description
* [reflected_signal.m](https://github.com/tuanna52/GPR-Corroded-Rebar-Detection/blob/main/reflected_signal.m): 
This file generates data to simulate the multiple reflections of EM wave from layer structure in two cases – pulse source and sinusoidal source.

* [foward_model_simulation.m](https://github.com/tuanna52/GPR-Corroded-Rebar-Detection/blob/main/foward_model_simulation.m):
This file contains the implementation of forward model and generates simulation data as well as Fourier transform.

* [Steel_In_Sandstone_Raw_Data_Only.mat](https://github.com/tuanna52/GPR-Corroded-Rebar-Detection/blob/main/Steel_In_Sandstone_Raw_Data_Only.mat):
This file contains the measured data from the conducted experiment.

* [data_processing.m](https://github.com/tuanna52/GPR-Corroded-Rebar-Detection/blob/main/data_processing.m):
This file pre-processes the data imported from “Steel_In_Sandstone_Raw_Data_Only.mat” file. Then applies inverse model to derive physical parameters of each layer inside a structure.

* [inverse_model.m](https://github.com/tuanna52/GPR-Corroded-Rebar-Detection/blob/main/inverse_model.m):
This file contains the function implementing the inverse model used in “data_processing.m” file.

* [fw_model.m](https://github.com/tuanna52/GPR-Corroded-Rebar-Detection/blob/main/fw_model.m):
This file contains the function implementing the forward model used in inverse_model() function as objective function for Levenberg-Marquardt curve-fitting algorithm.

* [simp_fw_model.m](https://github.com/tuanna52/GPR-Corroded-Rebar-Detection/blob/main/simp_fw_model.m):
This file contains the function implementing the simplified forward model used in inverse_model() function to optimize the permittivity and layers’ thickness before applying Levenberg-Marquardt algorithm.
