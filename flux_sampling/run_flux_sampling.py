import cobra
from cobra.sampling import sample
from cobra.sampling import OptGPSampler
import pandas as pd
import time
import glob
import os

model_path = '/home/maziya/CollateralLethality-2022new/models_xml_wo_outliers'
model_path = '\CollateralLethality-2022new\Breast'
os.chdir(model_path)
media = pd.read_csv('rpmi.csv', sep=',', header=0)
medium = media.set_index('metabolite').to_dict()
# medium = media.to_dict()
sampling_result = {}
sampling_result_valid ={}
model_files = glob.glob('*.xml')
for i in model_files:
    start = time.time()
    model = cobra.io.read_sbml_model(i)
    temp = medium.copy()
    for j in temp['flux'].copy():
        if j not in model.exchanges:
            temp['flux'].pop(j)
    model.medium.update(temp['flux'])
    model.exchanges.get_by_id('EX_glc(e)').upper_bound = -0.5
    optgp = OptGPSampler(model, processes=2)
    sampling_result[i] = [s for s in optgp.batch(100, 20)]
    sampling_result_valid[i] = []
    for sam in sampling_result[i]:
        s_valid = sam[optgp.validate(sam) == "v"]
        sampling_result_valid[i].append(len(s_valid))
    end = time.time()
    print('Time taken at the end of ', str(i), 'model: ', end-start)
final_dict = {}
for i in sampling_result:
    final_dict[i] = pd.concat(sampling_result[i])
    final_dict[i].T.to_csv(model_path+i.replace('.xml', '.csv'))
