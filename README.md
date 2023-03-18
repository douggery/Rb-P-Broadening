# Rb-P-Broadening
Pressure broadening fits of microcells fabricated using a new method

Silicon was micromachined with 100 cavities which were then encapsulated
using a new bonding technique using anodic bonding. 

The resulting vapor cells were measured for purity with broadening
of the spectroscopy being an indicator of impurities. 

In this code, the spectra of 100+ cells were measured and fit using
a deterministic model and removed systematic effects such as laser
inensity fluctuations, laser linewidth, and laser amplitude slope

Using the hyperfine features convolved with a voigt model, we reconstructed
the spectra correctly and were able to extract the features of the system
which were displayed in a histogram.

![WLD13_residual_pressures](https://user-images.githubusercontent.com/30641156/226086546-038cf344-53f3-4a38-a6d4-cf987f3605fa.png)

Earlier generation with 109 fits compared and aligned appropriately (laser frequency drift
was compesated in software) and my atrocious choice of colors was remedied.

![109 FITS with histogram zoom](https://user-images.githubusercontent.com/30641156/226086610-ec63e779-4792-4302-afe3-14bcec82990c.png)
