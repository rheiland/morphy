# morphy (for ws2022)- you know who you are


Latest:
```
make -j2
mv morph PhysiCell-model-builder-2.7.3   # want to have the executable in the GUI
cd PhysiCell-model-builder-2.7.3
python bin/pmb.py --studio    # or just try to run "studio.sh"
```
then in the PMB, select the Run tab, press Run Simulation. In Plot tab, Play.

Old/originally:
```
make -j2
morph
cd output
python cells_info.py 0
python anim_morph.py
```


Some edits I've made:
* fixed a bug in modules/PhysiCell_MultiCellDS.cpp
* "flattened" (filled in all params) of config/PhysiCell_settings.xml (but may have clobbered some of your desired params, e.g., I changed some of the original domain params)
* edited setup_tissue() in custom.cpp: for now, just define 4 cells; plus, I think the `axis_*` params were missing from <user_params>
* added `beta/anim_morph.py` for simple plotting (and also into /output, but it may get deleted if the `clean` or other Makefile target does so)
* cloned/copied the PMB into the root dir (might make it unique to this project)
* hacked `phenotype_function` in custom.cpp to do something silly, but at least provide dynamics, with the `axis_*` custom vars
* continue to make tweaks to anim_morph.py
* more to come...

Overnight:
* started using the PMB; added ellipse plotting to it
* the PMB's /data/morphy.xml is the latest config file
* added Didi's tweaks to setup_tissue() and new user_params to the .xml
* in vis_tab.py, ellipse's angle is a function of the cell's motility vector
