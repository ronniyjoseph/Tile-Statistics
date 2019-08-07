import os
from astropy.io import fits
import numpy
import subprocess
import multiprocessing
from functools import partial

def main(download = False):
    obsids_list = "Ultimate-EOR-obsids-2014-19-full.txt"
    metafits_folder = '/mnt/data/PPDs'
    local_folder = "/data/rjoseph/Beam_Perturbations/Dipole_Statistics/ppd_metafits_all_EoR/"
    if download:
        download_metafits_ppds(obsids_list, metafits_folder, local_folder)

    #count_broken_tiles_all_data(metafits_folder)
    count_broken_tiles_local_parallel(local_folder, file_name="All_EoR_Stats.txt", include_flagged=False)
    return


def count_broken_tiles_all_data(metafits_folder, include_flagged= False):
    print("Opening", metafits_folder)
    ppd_folder_list = os.listdir(metafits_folder)


    for ppd_folder in ppd_folder_list:
        print("Opening folder", ppd_folder)
        counter = 0
        obsid_list = os.listdir(metafits_folder + '/' + ppd_folder)

        tile_statistics = []
        print("Counting tiles with 1 broken dipole and with 2 broken dipoles")
        for file in obsid_list:
            filename_split = file.split("_")
            if filename_split[1] == "metafits":

                obsid_file_path = metafits_folder + "/" + ppd_folder + "/" + file

                hdu = fits.open(obsid_file_path)
                delay_table = hdu[1].data

                if include_flagged:
                    non_flagged_tile_indices = None
                else:
                    non_flagged_tile_indices = numpy.where(delay_table['Flag'] != 1)

                broken_tile_indices, broken_dipole_indices = numpy.where(delay_table['Delays'][non_flagged_tile_indices[0]] == 32)

                broken_tile_numbers = delay_table['Antenna'][non_flagged_tile_indices][broken_tile_indices]
                broken_tiles, tile_occurrence = numpy.unique(broken_tile_numbers, return_counts = True)

                single_pol_indices = tile_occurrence[tile_occurrence == 1]
                double_pol_indices = tile_occurrence[tile_occurrence == 2]

                metadata = [int(filename_split[0]), len(single_pol_indices), len(double_pol_indices)]
                tile_statistics.append(metadata)

                hdu.close()
                counter += 1

        numpy.savetxt("broken_tile_count_" + str(ppd_folder) + ".txt", numpy.array(tile_statistics), fmt='%i',
                      header="obsid 1dipole_count 2dipole_count")
    return

def count_broken_tiles_local(obsid_path, file_name = "Tile_Statistics.txt", include_flagged= False):
    print("Opening", obsid_path)
    obsid_list = os.listdir(obsid_path)
    tile_statistics = []
    counter = 0

    print("Counting tiles with 1 broken dipole and with 2 broken dipoles")
    for file in obsid_list:
        filename_split = file.split("_")

        obsid_metafits = obsid_path + "/" + file
        hdu = fits.open(obsid_metafits)
        if not counter % int(len(obsid_list)*0.1):
            print(f"{counter/len(obsid_list)*100}%")
        # the relevant data is in the 1 hdu table
        delay_table = hdu[1].data
        if include_flagged:
            non_flagged_tile_indices = None
        else:
            non_flagged_tile_indices = numpy.where(delay_table['Flag'] != 1)

            broken_tile_indices, broken_dipole_indices = numpy.where(delay_table['Delays'][non_flagged_tile_indices[0]] == 32)

            broken_tile_numbers = delay_table['Antenna'][non_flagged_tile_indices][broken_tile_indices]
            broken_tiles, tile_occurrence = numpy.unique(broken_tile_numbers, return_counts = True)

            single_pol_indices = tile_occurrence[tile_occurrence == 1]
            double_pol_indices = tile_occurrence[tile_occurrence == 2]

            if hdu[0].header['MODE'] == 'NO_CAPTURE':
                metadata = [int(filename_split[0]), numpy.nan, numpy.nan]
            else:
                metadata = [int(filename_split[0]), len(single_pol_indices), len(double_pol_indices)]
            tile_statistics.append(metadata)

            hdu.close()
            counter += 1
    numpy.savetxt(file_name, numpy.array(tile_statistics), fmt='%i', header="obsid 1dipole_count 2dipole_count")

    return

def count_broken_tiles_local_parallel(obsid_path, file_name = "Stats.txt", include_flagged= False):
    print("Opening", obsid_path)
    obsid_list = sorted(os.listdir(obsid_path))


    pool = multiprocessing.Pool(7)
    tile_statistics = pool.map(partial(single_count, obsid_path, obsid_list, include_flagged), range(len(obsid_list)))

    print("Counting tiles with 1 broken dipole and with 2 broken dipoles")

    numpy.savetxt(file_name, numpy.array(tile_statistics), fmt='%i', header="obsid 1dipole_count 2dipole_count")
    return


def single_count(obsid_path, obsid_list, include_flagged, index):


    filename_split = obsid_list[index].split("_")
    obsid_metafits = obsid_path + "/" + obsid_list[index]

    try:
        hdu = fits.open(obsid_metafits)
    except:
        metadata = [int(filename_split[0]), -1, -1]
    if not index % int(len(obsid_list) * 0.1):
        print(f"{index / len(obsid_list) * 100}%")

    # the relevant data is in the 1 hdu table
    delay_table = hdu[1].data

    if include_flagged:
        non_flagged_tile_indices = None
    else:
        non_flagged_tile_indices = numpy.where(delay_table['Flag'] != 1)

        broken_tile_indices, broken_dipole_indices = numpy.where(
            delay_table['Delays'][non_flagged_tile_indices[0]] == 32)

        broken_tile_numbers = delay_table['Antenna'][non_flagged_tile_indices][broken_tile_indices]
        broken_tiles, tile_occurrence = numpy.unique(broken_tile_numbers, return_counts=True)

        single_pol_indices = tile_occurrence[tile_occurrence == 1]
        double_pol_indices = tile_occurrence[tile_occurrence == 2]

        if hdu[0].header['MODE'] == 'NO_CAPTURE':
            metadata = [int(filename_split[0]), numpy.nan, numpy.nan]
        else:
            metadata = [int(filename_split[0]), len(single_pol_indices), len(double_pol_indices)]

        hdu.close()
        return metadata

def count_broken_tiles(obsids_list, metafits_folder, include_flagged= False):
    print("Opening", obsids_list)
    obsids_array = numpy.loadtxt(obsids_list, ndmin = 1)

    if len(obsids_array) > 1:
        obsids_array = numpy.sort(obsids_array)

    counter = 0
    tile_statistics = numpy.zeros((len(obsids_array), 3), dtype = int)
    dipole_statistics = numpy.zeros((len(obsids_array), 1 + 16), dtype = int)

    tile_count1_file = open("broken_tile_ids_1dipole_" + obsids_list, "w")
    tile_count2_file = open("broken_tile_ids_2dipole_" + obsids_list, "w")

    print("Counting tiles with 1 broken dipole and with 2 broken dipoles")
    for obsid in obsids_array:
        obsid_path = str(int(obsid))[:4]
        obsid_metafits = metafits_folder + "/" + obsid_path + "/" + str(int(obsid)) + "_metafits_ppds.fits"
        hdu = fits.open(obsid_metafits)

        #the relevant data is in the 1 hdu table
        delay_table = hdu[1].data
        if include_flagged:
            non_flagged_tile_indices = None
        else:
            non_flagged_tile_indices = numpy.where(delay_table['Flag'] != 1)

        broken_tile_indices, broken_dipole_indices = numpy.where(delay_table['Delays'][non_flagged_tile_indices[0]] == 32)

        broken_tile_numbers = delay_table['Antenna'][non_flagged_tile_indices][broken_tile_indices]
        broken_tiles, tile_occurrence = numpy.unique(broken_tile_numbers, return_counts = True)

        single_pol_indices = tile_occurrence[tile_occurrence == 1]
        double_pol_indices = tile_occurrence[tile_occurrence == 2]

        tile_statistics[counter, 0] = obsid
        tile_statistics[counter, 1] = len(single_pol_indices)
        tile_statistics[counter, 2] = len(double_pol_indices)

        # Extract all broken tile names per obsid
        all_tile_names_sorted = delay_table['TileName'][::2][delay_table['Antenna'][::2].argsort()]
        tile_names_1dipole = all_tile_names_sorted[broken_tiles[tile_occurrence == 1]]
        tile_names_2dipole = all_tile_names_sorted[broken_tiles[tile_occurrence == 2]]

        tile_count1_file.write(str(int(obsid)) + " " + ' '.join([name for name in tile_names_1dipole]) + "\n")
        tile_count2_file.write(str(int(obsid)) + " " + ' '.join([name for name in tile_names_2dipole]) + "\n")


        # Extract dipole information
        broken_dipole_numbers = numpy.arange(0,16,1)[broken_dipole_indices]
        broken_dipoles, dipole_occurence = numpy.unique(broken_dipole_numbers, return_counts = True)

        dipole_statistics[counter, broken_dipoles + 1] = dipole_occurence

        hdu.close()
        counter += 1

    tile_count1_file.close()
    tile_count2_file.close()


    numpy.savetxt("broken_tile_count_" + obsids_list, tile_statistics, fmt ='%i', header  = "obsid 1dipole_count 2dipole_count" )
    numpy.savetxt("broken_dipole_count_"+ obsids_list, numpy.sum(dipole_statistics, axis = 0)/len(obsids_array))
    return


def download_metafits_ppds(obsids_list, metafits_folder, local_folder):
    print("Downloading Files from", obsids_list)
    obsids_array = numpy.loadtxt(obsids_list, ndmin = 1)

    if len(obsids_array) > 1:
        obsids_array = numpy.sort(obsids_array)

    for obsid in obsids_array:
        obsid_path = str(int(obsid))[:4]
        obsid_metafits = metafits_folder + "/" + obsid_path + "/" + str(int(obsid)) + "_metafits_ppds.fits"

        process = subprocess.Popen(["scp", obsid_metafits, "Danakat:" + local_folder])
        process.communicate()
    return









if __name__ == "__main__":
    main()