# Realtime-Diffraction-PBR

This plugin was developed for the publication "A Practical Real-Time Model for Diffraction on Rough Surfaces" - Clausen et al. 2024
link-to-publication

Insert Teaser Image

Description of repository. Developed and tested with Unity 2021.3.6f1, Universal RP 12.1.7, OpenXR Plugin 1.4.2

# Ultra Quick Start
- Download repository and open as Unity project (tested on Unity 2021.3.6f1)
- Import test scene package from Resources link
- Open the test scene
- For VR to work properly in the test scene, the "Starter Assets" from the xr.interaction.toolkit package have to be imported via the Unity package manager

# Integration into existing project
- The diffraction material implementation is based on the Universal Rendering Pipeline 12.1.7. We copied the minimum set of modified shaders into this repository. The rest of the shaders are assumed to come from your local URP 12.1.7 package. It may work on other URP versions as well, however we have not tested it, or validated the visual result
- Simply assign the "CustomShaders/DiffractionLit" material to any object with an existing UV mapping

# Diffraction Material Settings
- Todo

# Known Limitations
- Todo

# Possible Future Improvements
- Todo
- 
# Resources
Museum test scene with some materials
https://1drv.ms/u/s!Ap1NX8WBfJHQg6wuDvZ5oTYSAKzumw
