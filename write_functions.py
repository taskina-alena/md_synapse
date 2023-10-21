def write_config(configName, system):
    if system.bi_domain:
        n_bond_types = system.num_v + 1
        n_angles_types = 1
    else:
        n_bond_types = system.num_v
        n_angles_types = 0
    header = ["LAMMPS Description \n \n",
              "\t " + str(system.numAll) + " atoms \n \t " + str(system.num_bonds) + " bonds \n \t " +
              str(system.num_angles) + " angles \n \t 0 dihedrals \n \t 0 impropers \n",
              "\n \t " + str(system.num_types) + ' atom types \n \t ' + str(n_bond_types) +
              ' bond types \n \t ' + str(n_angles_types) + ' angle types \n \t 0 dihedral types \n \t 0 improper types \n',
              "\n \t " + str(-system.Lx * 0.5) + " " + str(system.Lx * 0.5) + " xlo xhi\n \t",
              str(-system.Ly * 0.5) + " " + str(system.Ly * 0.5) + " ylo yhi \n \t",
              str(-system.Lz * 0.5) + " " + str(system.Lz * 0.5) + " zlo zhi\n"]

    header.append("\nMasses \n \n")

    for i in range(len(system.type_mass_list)):
        header.append("\t %i %.4f \n" % (system.type_mass_list[i][0], system.type_mass_list[i][1]))

    f = open(configName, "w")

    for item in header:
        f.write(item)

    for item in system.coords:
        f.write(item)

    for item in system.bonds:
        f.write(item)

    if system.bi_domain:
        for item in system.angles:
            f.write(item)

    f.close()


def write_input(input_file, config_file, system, folderPattern, filePattern, seed, dump_step, run_steps):

    f = open(input_file, "w")

    f.write("units                  lj                                      \n")
    f.write("dimension 2\n")
    f.write("atom_style             full  \n")
    f.write("boundary               p p p   \n")
    f.write(f"read_data {config_file} \n")

    f.write("neighbor               0.3 bin\n")
    f.write("neigh_modify           every 1 delay 1\n")

    for item in system.bonds_interactions:
        f.write(item)

    if system.bi_domain:
        for item in system.angles_interactions:
            f.write(item)

    f.write(
        '''
        ##======================================##
        ## Interactions
        ##======================================##
        \n''')

    for item in system.interactions:
        f.write(item)

    f.write(
        '''
        ##======================================##
        ## Setting up groups
        ##======================================##
        \n''')
    # f.write("group                  head type 1\n")

    f.write(
        '''
        ##======================================##
        ## Computes
        ##======================================##
        \n''')

    f.write(
        '''
        ##======================================##
        ## Fixes
        ##======================================##
        \n''')

    # f.write('minimize 1.0e-4 1.0e-6 1000 10000 \n')

    f.write("velocity               all create 1.0 1 \n")

    f.write("log                    %s/Log_%s.dat\n" % (folderPattern, filePattern))

    f.write("fix                    fLANG all langevin 1.0 1.0 1.0 %i\n" % seed)
    f.write("fix                    fNVE all nve\n")
    f.write(
        '''
        ##======================================##
        ## Analysis
        ##======================================##
        \n''')

    f.write("thermo                 %i\n" % dump_step)
    f.write("thermo_style           custom step temp press etotal epair\n")

    f.write("thermo_modify          flush yes\n")
    f.write("timestep               0.005\n")

    f.write(
        '''
        ##======================================##
        ## Equilibration before output
        ##======================================##
        \n''')

    f.write(
        '''
        ##======================================##
        ## Run with output
        ##======================================##
        \n''')
    f.write("fix enforce_2d all enforce2d\n")
    f.write("dump                       1 all custom %i %s/Movie_%s.xyz id type mol x y z q  \n" % (
        dump_step, folderPattern, filePattern))

    f.write("run                    %i\n" % run_steps)
    f.close()