### Optical Modeling of Superconducting Nanowire Single Photon Detectors

This repository contains MATLAB scripts that I wrote to model the absorptance of light in superconducting nanowire single photon detectors (SNSPDs) as well as the reflectance, transmittance and total absorptance of a structure. `all_optical_simulations.m` has the basic script for simulating illumination through a substrate at normal incidence, and the other scripts listed below modify it to explore other aspects of the device geometry or physics.

* `all_optical_simulations_varying_wavelength.m` gives the optical characteristics of an SNSPD at different wavelengths and relies on external files with the optical constants.
* `all_optical_simulations_coherent_substrate.m` considers the case where the substrate thickness is *not* greater than the coherence length of light (e.g., when a CW laser is used).
* `all_optical_variable_angle.m` varies the angle of incidence and was used to model SNSPDs on SiNx membranes on cleaved optical fibers.
* `all_optical_simulations_2D_maps.m` produces two-dimensional "heat maps" that show the absorptance as the thickness and fill factor of the devices are varied.
