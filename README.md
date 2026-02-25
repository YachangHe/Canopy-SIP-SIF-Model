# Canopy-SIP SIF Model

contains the MATLAB implementation of the canopy-SIP SIF model



\# Canopy-SIP SIF Model



This repository contains the MATLAB implementation of the canopy-SIP SIF model, which simulates canopy solar-induced chlorophyll fluorescence (SIF) anisotropy. 



The model was developed by combining geometric-optical (GO) theory to account for asymmetric leaf SIF forward and backward emissions at the first-order scattering, and by modeling multiple scattering based on the spectral invariants theory (p-theory). This approach avoids the dependence on complex 3D radiative transfer models while maintaining high accuracy over 3D canopy structures.



\## Repository Structure



\* `demo.m`: The main script to demonstrate the simulation of SIF anisotropy and generate the output spectrum.

\* `/CI\_2/`: Contains structural data, including total gap fraction, within-crown gap fraction, and clumping index.

\* `leafopt\_P6.mat`, `Esun\_SIP.mat`, `rs\_01.mat`: Input spectral data (leaf optics, solar irradiance, and soil reflectance).

\* `HET01\_true\_all.txt`: Scene information used for hemispherical interceptance calculation.

\* `\*.m`: Core sub-functions for LIDF generation, scattering probabilities calculation, and radiative transfer simulation.



\## System Requirements



\* \*\*Software\*\*: MATLAB (Tested on R2023a or later. Previous versions may work but are not strictly tested).

\* \*\*Toolboxes\*\*: No specific external toolboxes are required. The code relies on standard MATLAB functions.



\## Usage



1\.  Clone or download this repository to your local machine.

2\.  Open MATLAB and set the repository folder as your Current Folder.

3\.  Open and run the `demo.m` script.

4\.  The script will automatically load the necessary structural and optical data, perform the SIF directional emission calculations, save the results to `demo\_result\_2026.mat`, and plot the simulated Canopy SIF Spectrum at the nadir view.



\## Citation



If you use this model or code in your research, please cite the following paper:



> He, Yachang, Yelu Zeng, Dalei Hao, Nikolay V. Shabanov, Jianxi Huang, Gaofei Yin, Khelvi Biriukova et al. "Combining geometric-optical and spectral invariants theories for modeling canopy fluorescence anisotropy." Remote Sensing of Environment 323 (2025): 114716. https://doi.org/10.1016/j.rse.2025.114716


\## Model Heritage & Acknowledgments

This model represents a continuous effort in canopy radiative transfer modeling. We gratefully acknowledge the developers of the following models, whose foundational work, theoretical frameworks, and validation tools significantly contributed to this project:

* **Original SIP Model**: This model is built upon the foundational work of the FluorRTER model model.
    > Zeng, Yelu, Grayson Badgley, Min Chen, Jing Li, Leander DL Anderegg, Ari Kornfeld, Qinhuo Liu et al. "A radiative transfer model for solar induced fluorescence using spectral invariants theory." Remote Sensing of Environment 240 (2020): 111678.
* **PATH_RT Model**: Parts of the sub-functions used for calculating the four-component gap fractions and hotspot effects are derived from the PATH_RT model developed by Dr. Weihua Li and colleagues.
    > Li, W., Yan, G., Mu, X., Tong, Y., Zhou, K., & Xie, D. (2024). Modeling the hotspot effect for vegetation canopies based on path length distribution. *Remote Sensing of Environment*, 303, 113985.
* **LESS Model**: The comparison of our discrete canopy model simulations was conducted using the LESS framework.
    > Qi, J., Xie, D., Yin, T., Yan, G., Gastellu-Etchegorry, J. P., Li, L., Zhang, W., Mu, X., & Norford, L. K. (2019). LESS: LargE-Scale remote sensing data and image simulation framework over heterogeneous 3D scenes. *Remote Sensing of Environment*, 221, 695-706.

\## License

This project is licensed under the MIT License - see the LICENSE file for details.

