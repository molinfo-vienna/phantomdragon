import MDAnalysis as mda
import MDAnalysis.transformations as trans

u = mda.Universe("/data/shared/projects/master_linda/1bcd/openmm/step3_input.psf",["/data/shared/projects/master_linda/1bcd/openmm/step5_1.dcd","/data/shared/projects/master_linda/1bcd/openmm/step5_2.dcd","/data/shared/projects/master_linda/1bcd/openmm/step5_3.dcd","/data/shared/projects/master_linda/1bcd/openmm/step5_4.dcd","/data/shared/projects/master_linda/1bcd/openmm/step5_5.dcd","/data/shared/projects/master_linda/1bcd/openmm/step5_6.dcd","/data/shared/projects/master_linda/1bcd/openmm/step5_7.dcd","/data/shared/projects/master_linda/1bcd/openmm/step5_8.dcd","/data/shared/projects/master_linda/1bcd/openmm/step5_9.dcd","/data/shared/projects/master_linda/1bcd/openmm/step5_10.dcd"],in_memory=True)

print("loading complete")

protein = u.select_atoms('protein')
not_protein = u.select_atoms('not protein')


transforms = [trans.unwrap(u.atoms),
              trans.center_in_box(protein, center='geometry'),
              trans.wrap(not_protein,compound='residues')]

u.trajectory.add_transformations(*transforms)

print("transform complete")
    
all = u.select_atoms('all')

with mda.Writer("test.dcd", all.n_atoms) as W:
    for ts in u.trajectory:
        W.write(all)