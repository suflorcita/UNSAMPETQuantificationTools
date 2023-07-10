import os
def write_line_quantification(subject, RMN_path, FDG_path, output_path, path_file):
    with open(path_file, 'a') as f:
        f.write(
            f'python3 main.py '
            f'{FDG_path} '
            f'{RMN_path} '
            f'{subject} '
            f'{output_path}\n'
        )
    pass

def write_line_chart(subject, output_chart, subject_folder, path_file):
    with open(path_file, 'a') as f:
        f.write(
            f'python3 charts.py '
            f'{output_chart} '
            f'{subject_folder} '
            f'{subject} \n'
        )
    pass

def generate_script_pet_quant(output_directory, script_path):
    files = sorted(os.listdir(output_directory))

    # clear the data in the FreeSurfer file and write first line
    with open(script_path, 'w') as file:
        file.write('#!/bin/bash\n')

    for subject in files:
        dir_subject = output_directory + "/" + subject


        # Nifti path
        nifti_subject = dir_subject + "/Nifti"

        if os.path.isdir(nifti_subject):
            # RMN path
            try:
                RMN_path = nifti_subject + "/RMN/" + os.listdir(nifti_subject + "/RMN")[0]
            except IndexError:
                print(f"{subject} ...")
                continue
            # PET path
            PET_path = nifti_subject + "/FDG-PET"

            write_line_quantification(subject, RMN_path, PET_path,output_directory, script_path)

def generate_script_charts(output_directory, script_path):
    files = sorted(os.listdir(output_directory))

    # clear the data in the FreeSurfer file and write first line
    with open(script_path, 'w') as file:
        file.write('#!/bin/bash\n')

    for subject in files:

        if subject[0] != '0':
            print(subject)
            continue
        dir_subject = output_directory + "/" + subject

        # Chart path
        chart_path = dir_subject + "/charts"

        write_line_chart(subject, chart_path, dir_subject, script_path)



if __name__ == '__main__':

    output_directory = "/home/sol/PET_MRI/Procesamiento/Procesar"
    output_directory_chart = "/home/sol/PET_MRI/Procesamiento"

    # create a script to run FreeSurfer
    script_path_quant = "/home/sol/PET_MRI/quantification.sh"
    script_path_chart = "/home/sol/PET_MRI/charts.sh"

    generate_script_pet_quant(output_directory, script_path_quant)
    #generate_script_charts(output_directory_chart, script_path_chart)




    #
