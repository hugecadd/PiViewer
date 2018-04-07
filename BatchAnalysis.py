from PiViewer import find_PiPi
import os


if __name__ == '__main__':

    if False:   # test if single pdb works
        # pdb_file = r'D:\Iridium_HT_deposited\1a28_d1refined.pdb'
        pdb_file = r'C:\CloudStation\Epicat\Git\PiViewer\1ACJ.pdb'
        lig_name = 'THA'
        print find_PiPi(pdb_file, lig_name, 5.0, 25, 80, 0)

    # Parse the list file containing PDB codes and ligand res names
    pdb_list_file = r'.\Dataset\Iridium_HT_PDB_list.txt'
    with open(pdb_list_file, 'r') as fin:
        os.chdir('C:\CloudStation\Epicat\Git\PiViewer\Dataset\Iridium_HT_deposited')

        for line in fin:
            items = line.split()
            pdb_code = items[0]
            lig_list = items[1].split(',')
            print pdb_code,
            # print lig_list,

            # generate pdb file path
            pdb_file = pdb_code.lower() + '_d1refined.pdb'

            # in many cases the ligand res name is UNL instead of specified in the list
            lig_list.append('UNL')

            # detect Pi-Pi interactions
            total_found = 0
            total_unfound = 0
            for lig_name in lig_list:
                count = find_PiPi(pdb_file, lig_name, verbose=0)
                if count > 0:
                    total_found += count
                elif count == -1:
                    total_unfound += 1
            if total_unfound == len(lig_list):
                print -1
            else:
                print total_found
    #