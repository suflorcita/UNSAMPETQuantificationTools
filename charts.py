import matplotlib.pyplot as plt
import pandas as pd
import numpy as np 
import SimpleITK as sitk
import sys
from scipy import ndimage as ndi



def mean_dataframe(df_intensity_region):
    '''Creates a dataframe with the intensity of the signal for region in the atlas.
            Args:
               df_intensity_region: 
               A dataframe with the intensity of the signal for region in each hemisphere.

            Returns:
                df_mean_regions: 
                A dataframe with the mean intensity of the signal for region in both hemispheres.
            '''

    # mean of the two regions
    df_mean_regions = pd.DataFrame(df_intensity_region.groupby(["structure"])["normalization"].mean())

    # error
    errors = []

    # Iterate over the dataframe
    for index, row in df_mean_regions.iterrows():

        # encuentro las regiones que coinciden con el índice
        structures = df_intensity_region.loc[df_intensity_region["structure"] == index]

        if len(structures) != 2:
            error = np.nan
            errors.append(error)
            continue

        structure_left = structures.loc[structures["hemisphere"] == "L"]
        structure_right = structures.loc[structures["hemisphere"] == "R"]

        value_left = float(structure_left["normalization"])
        value_right = float(structure_right["normalization"])
        value_mean = row["normalization"]

        error = abs(((value_left - value_right) / value_mean) * 100)
        errors.append(error)

    df_mean_regions["error"] = errors

    df_mean_regions = df_mean_regions.dropna()

    return df_mean_regions


def plot_images(T1, cambios_cerebro, n): 
    img_T1 = sitk.GetArrayFromImage(T1)
    img_cambios_cerebro = sitk.GetArrayFromImage(cambios_cerebro)


    fig_rows = 4
    fig_cols = 4
    n_subplots = fig_rows * fig_cols
    n_slice = img_cambios_cerebro.shape[0]
    step_size = n_slice // n_subplots
    plot_range = n_subplots * step_size
    start_stop = int((n_slice - plot_range) / 2)

    fig, axs = plt.subplots(fig_rows, fig_cols)

    for idx, img in enumerate(range(start_stop, plot_range, step_size)):
        axs.flat[idx].imshow(img_T1[img, :, :], cmap='gray', alpha=1)
        
        if n == 0: 
            axs.flat[idx].imshow(img_cambios_cerebro[img, :, :], cmap='bwr', vmin=-50, vmax=50, alpha=0.5)
        elif n == 1: 
            axs.flat[idx].imshow(ndi.rotate(img_cambios_cerebro[:, img, :], 180), cmap='bwr', vmin=-50, vmax=50, alpha=0.5)
        elif n == 2: 
            axs.flat[idx].imshow(ndi.rotate(img_cambios_cerebro[:, :, img], 180), cmap='bwr', vmin=-50, vmax=50, alpha=0.5)
        axs.flat[idx].axis('off')

    

def intensity_regions_bar_chart(regions, intensity1, intensity2, colors=None, hemisphere=None, title=None):
    """
        Creates a bar chart using the brain regions in the x-axis and the signal intensity in the y-axis.

        Parameters
        ----------
        regions : list of brain regions
        intensity: list of signal intensity for every region in regions.
        title: title of the graph

    """

    regions_graph = [str(i) + "_" + regions for i, regions in enumerate(regions)]
    if hemisphere:
        regions_graph= [region + "_left" if i % 2 == 0 else region + "_right"  for i, region in enumerate(regions_graph)]
    y_pos = np.arange(len(regions_graph))
    intensity_pet = [intensity1, intensity2]

    plt.bar(y_pos-0.2, intensity_pet[0], color='orange', width= 0.4, label=regions_graph)
    plt.bar(y_pos+0.2, intensity_pet[1], color='cornflowerblue', width= 0.4)

    plt.xlabel("Regiones", size=16)
    plt.ylabel("Intensidad", size=16)
    plt.xticks([i for i in range(len(regions))])

    if title != None: plt.title(title, size=16)
    else:
        plt.title("Gráfico por intensidades", size=16)
    plt.legend(bbox_to_anchor=(1.001, 1), fontsize='small')
    plt.subplots_adjust(right=0.80)



if __name__ == '__main__':

    # images
    if len(sys.argv) > 1:
        output_charts = sys.argv[1]
        CSV_folder = sys.argv[2]
        subject_name = sys.argv[3]
    else:
        output_charts = "../Procesado/SubjectCEUNIM/charts"
        CSV_folder = "../Procesado/SubjectCEUNIM/CSV"
        subject_name = "SubjectCEUNIM"
    # Bar chart mean for brain structure


    # Read CSV MNI152 and subject into a dataframe 
    MNI152 = pd.read_csv(CSV_folder + "/" + "MNI152_intensity.csv")
    subject = pd.read_csv(CSV_folder + "/" + "subject_intensity.csv")

    init = 0
    step = len(MNI152) // 3 + 1
    finish = len(MNI152)
    j = 1

    # Regions
    for i in range(init, finish, step):
        regions = MNI152[i:i + step]["structure"]
        hemispheres = MNI152[i: i + step]["hemisphere"]
        name_regions = pd.Series(regions + "_" + hemispheres)

        # intensity values for MNI152 and subject
        intensity_MNI152 = MNI152[i: i + step]["normalization"]
        intensity_subject = subject[i: i + step]["normalization"]

        chart_df = pd.DataFrame()
        chart_df["name_regions"] = name_regions
        chart_df["intensity_MNI152"] = intensity_MNI152
        chart_df["intensity_subject"] = intensity_subject
        chart_df = chart_df.dropna()

        name_chart = f"Intensidad por región normalizado al cerebelo {j}"

        plt.figure(j)
        intensity_regions_bar_chart(chart_df["name_regions"], chart_df["intensity_subject"], chart_df["intensity_MNI152"],
                                    title=name_chart)
        name_chart = f'{output_charts}/{subject_name}_all_regions_{j}.png'

        # Set figure size
        fig = plt.gcf()  # Get the current figure
        fig.set_size_inches(16, 8)  # Set the size in inches (width, height)
        plt.savefig(name_chart, dpi=600)

        j += 1

    # Mean for region 
    df_MNI152_mean = mean_dataframe(MNI152)
    df_subject_mean = mean_dataframe(subject)

    regions_mean = df_MNI152_mean.index
    
    # Intensity values
    intensity_mean_Subject_normalization = df_subject_mean["normalization"]
    intensity_mean_MNI152_normalization = df_MNI152_mean["normalization"]
    
    # Error: 
    intensity_mean_Subject_error = df_subject_mean["error"]
    intensity_mean_MNI152_error = df_MNI152_mean["error"]

    j += 1
    plt.figure(j)
    intensity_regions_bar_chart(regions_mean, intensity_mean_Subject_normalization, intensity_mean_MNI152_normalization,
                                title="Grafico por intensidades, promedio por región")
    name_chart = f'{output_charts}/{subject_name}_mean_intensity.png'

    # Set figure size
    fig = plt.gcf()  # Get the current figure
    fig.set_size_inches(16, 8)  # Set the size in inches (width, height)
    plt.savefig(name_chart, dpi=600)

    j += 1
    plt.figure(j)
    intensity_regions_bar_chart(regions_mean, intensity_mean_Subject_error, intensity_mean_MNI152_error,
                                title="%error")
    name_chart = f'{output_charts}/{subject_name}_mean_error.png'

    # Set figure size
    fig = plt.gcf()  # Get the current figure
    fig.set_size_inches(16, 8)  # Set the size in inches (width, height)
    plt.savefig(name_chart, dpi=600)

    # Charts: Top ren regions more hypometabolic  
    
    # Read CSV: differences to MNI152
    df_diff = pd.read_csv(CSV_folder + "/" + "Intensity_image_changes.csv")
    
    # Sort the dataframe for change column in descendent order 
    df_diff_sorted = df_diff.sort_values("cambio")

    # Top Ten values in the Dataframe 
    top_ten = df_diff_sorted.head(10)

    # List of labels of top ten regions
    labels_top_ten = list(top_ten["n_label"])

    # Create a dataframe with top ten regions and the intensity value for MNI152
    top_ten_MNI152 = pd.DataFrame()

    for label in labels_top_ten:
        row = MNI152.loc[MNI152['n_label'] == label]
        top_ten_MNI152 = pd.concat([top_ten_MNI152, row], ignore_index=True)

    # Create a dataframe with top ten regions and the intensity value for Subject
    top_ten_Subject = pd.DataFrame()

    for label in labels_top_ten:
        row = subject.loc[subject['n_label'] == label]
        top_ten_Subject = pd.concat([top_ten_Subject, row], ignore_index=True)

    # Intensity values
    intensity_top_ten_subject_normalization = top_ten_Subject["normalization"]
    intensity_top_ten_MNI152_normalization = top_ten_MNI152["normalization"]

    # Regions
    regions_top_ten = top_ten_Subject["structure"]

    j += 1
    plt.figure(j)
    intensity_regions_bar_chart(regions_top_ten, intensity_top_ten_subject_normalization,
                                intensity_top_ten_MNI152_normalization,
                                title="Top 10 regions")
    name_chart = f'{output_charts}/{subject_name}_top_ten_regions.png'

    # Set figure size
    fig = plt.gcf()  # Get the current figure
    fig.set_size_inches(16, 8)  # Set the size in inches (width, height)
    plt.savefig(name_chart, dpi=600)

    print(f"charts {subject_name}:ok")



    # # Charts
    #
    # regions = regions_10_MNI_152["structure"]
    # hemispheres = regions_10_MNI_152["hemisphere"]

    # name_regions = []

    # for i, region in enumerate(regions):
    #     name_region = region + "-" + hemispheres[i]
    #     name_regions.append(name_region)

    # intensity_10_Subject = regions_10_MNI_152["normalization"]
    # intensity_10_MNI152 = regions_10_Subject["normalization"]

    # plt.savefig(output_graphs + "/graph2.png")

    # plt.figure(3)

    # intensity_regions_bar_chart(name_regions, intensity_10_Subject, intensity_10_MNI152,
    #                             title="Regions de mayor cambio")

    # plt.savefig(output_graphs + "/graph3.png")

    # # cambios cerebro
    # img_cambios_cerebro = sitk.GetArrayFromImage(cambios_cerebro)

    # # read T1 normalizada
    # path_T1 = "/home/sol/PET_MRI/Subject/Procesado/PET/Subject_to_MNI_152_ONLY_BRAIN/T1_MNI152_1mm_onlybrain_ANTsWarped.nii.gz"
    # sitk_T1 = sitk.ReadImage(path_T1)  # Read MNI 152
    # sitk_T1 = resample_sitk(sitk_T1, sitk_atlas_hammers)
    # img_T1 = sitk.GetArrayFromImage(sitk_T1)

    # sitk.WriteImage(sitk_T1, output_subject_PET + "/T1_resampleado_only_brain.nii.gz")

    # fig_rows = 4
    # fig_cols = 4
    # n_subplots = fig_rows * fig_cols
    # n_slice = img_cambios_cerebro.shape[0]
    # step_size = n_slice // n_subplots
    # plot_range = n_subplots * step_size
    # start_stop = int((n_slice - plot_range) / 2)

    # fig, axs = plt.subplots(fig_rows, fig_cols)

    # for idx, img in enumerate(range(start_stop, plot_range, step_size)):
    #     axs.flat[idx].imshow(img_T1[img, :, :], cmap='gray', alpha=1)
    #     axs.flat[idx].imshow(img_cambios_cerebro[img, :, :], cmap='bwr', vmin=-50, vmax=50, alpha=0.5)
    #     axs.flat[idx].axis('off')

    # plt.savefig(output_graphs + "/graph4.png")

    # plt.figure(4)

    # n_slice = img_cambios_cerebro.shape[1]
    # step_size = n_slice // n_subplots
    # plot_range = n_subplots * step_size
    # start_stop = int((n_slice - plot_range) / 2)
    # fig, axs = plt.subplots(fig_rows, fig_cols)

    # for idx, img in enumerate(range(start_stop, plot_range, step_size)):
    #     axs.flat[idx].imshow(ndi.rotate(img_T1[:, img, :], 180), cmap='gray', alpha=1)
    #     axs.flat[idx].imshow(ndi.rotate(img_cambios_cerebro[:, img, :], 180), cmap='bwr', vmin=-50, vmax=50, alpha=0.5)
    #     axs.flat[idx].axis('off')

    # plt.savefig(output_graphs + "/graph5.png")

    # plt.figure(5)

    # n_slice = img_cambios_cerebro.shape[2]
    # step_size = n_slice // n_subplots
    # plot_range = n_subplots * step_size
    # start_stop = int((n_slice - plot_range) / 2)
    # fig, axs = plt.subplots(fig_rows, fig_cols)

    # for idx, img in enumerate(range(start_stop, plot_range, step_size)):
    #     axs.flat[idx].imshow(ndi.rotate(img_T1[:, :, img], 180), cmap='gray', alpha=1)
    #     axs.flat[idx].imshow(ndi.rotate(img_cambios_cerebro[:, :, img], 180), cmap='bwr', vmin=-100, vmax=100, alpha=0.5)
    #     axs.flat[idx].axis('off')

    # plt.savefig(output_graphs + "/graph6.png")
    # plt.figure(6)

    # plt.tight_layout()
    # # plt.show()

