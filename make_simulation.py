from linker_class import LinkedVesicles
import os
import numpy as np
import write_functions


def create_files_structure(eps, num_vesicles, syn_density, type):
    filePattern = f"2d_{type}_eps{eps}_nv{num_vesicles}_sd{syn_density}"
    folderPattern = f"Results"
    gen_folders = ['Input/Scripts', 'Input/Configuration', folderPattern]

    for folder in gen_folders:
        if not os.path.exists(folder):
            os.makedirs(folder)

    return folderPattern, filePattern

if __name__ == "__main__":

    #eps_range = np.linspace(1, 6, 10)
    #lambda_range = np.linspace(0, 0.7, 10)

    eps_range = [6]
    lambda_range = [0.2]

    for bi_domain in [True, False]:
        for eps in eps_range:
            for syn_density in lambda_range:

                num_vesicles = 200
                sim_seed = 5
                real = 0
                dump_step = 1e4
                run_steps = 1e7
                side_box = 400

                lj_factor = 2 ** (1 / 6)
                global_cutoff = 3.0
                attr_range = 1

                kBond_vesicle = 100.0
                kBond_itra_syn = 100.0
                kBend = 1

                radii_distr = [4.17, 1.08]
                radii_vesicles = np.abs(np.random.normal(radii_distr[0], radii_distr[1], num_vesicles))
                radii = {'vesicles': radii_vesicles, 'synapsin': 0.5}

                def calc_syn_num(r, syn_density = syn_density):
                    cf = 2 * np.pi * r
                    return int(cf*syn_density)

                nums_syn = [calc_syn_num(r) for r in radii_vesicles]

                system = LinkedVesicles(num_vesicles, nums_syn, side_box, radii, bi_domain, kBend, kBond_vesicle, kBond_itra_syn, eps)

                if bi_domain:
                    folderPattern, filePattern = create_files_structure(eps, num_vesicles, syn_density, 'bi')
                else:
                    folderPattern, filePattern = create_files_structure(eps, num_vesicles, syn_density, 'mono')
                configName = f"Input/Configuration/Config_{filePattern}.dat"
                inputName = f"Input/Scripts/Input_{filePattern}.in"

                write_functions.write_config(configName, system)
                write_functions.write_input(inputName, configName, system, folderPattern, filePattern, sim_seed, dump_step, run_steps)





