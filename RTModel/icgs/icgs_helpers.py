import numpy as np
import pandas as pd
import glob

# modelcodes = ['PS', 'PX', 'BS', 'BO', 'LS', 'LX', 'LO', 'LK', 'TS', 'TX']
class FindFixedParameters:
    def __init__(self,eventname,modnumber,nfilt,grid_dictionary):
        self.eventname = eventname
        self.modnumber = modnumber
        self.nfilt = nfilt
        self.modelcodes =  ['PS', 'PX', 'BS', 'BO', 'LS', 'LX', 'LO', 'LK', 'TS', 'TX','TO']
        self.model_code = self.modelcodes[self.modnumber]
        self.grid_dictionary = grid_dictionary

        predecessor_model_map = {'LS':'PS',
                                 'LX':'PX',
        }

        if self.model_code not in predecessor_model_map.keys:
            raise ValueError(f"No predecessor model defined for {self.model_code}")
        run_map = {
        'PS': self.PS_InitConds,
        'PX': self.PX_InitConds,
        'BS': self.BS_InitConds,
        'BO': self.BO_InitConds,
        'LS': self.LS_InitConds,
        'LX': self.LX_InitConds,
        'LO': self.LO_InitConds,
        'LK': self.LK_InitConds,
        'TS': self.TS_InitConds,
        'TX': self.TX_InitConds,
        'TO': self.TO_InitConds
        }

        run_map[predecessor_model_map[self.model_code]]()

    def get_PS(self):
        model_list = glob.glob(f"{self.eventname}/Models/PS*")
        model_parameters = []
        for i in range(len(model_list)):
            modelfile = model_list[i]
            with open(modelfile) as f:
                line = f.readline()
                line = np.array(map(float,line.split(' ')))[[0,1,2,3,2*self.nfilt+4]]
                model_parameters.append(line)

        model_parameters = pd.DataFrame(model_parameters, columns = ['u0','tE','t0','rho','chi2'])
        model_parameters = model_parameters.sortby('chi2')


    def get_PX(self):
        raise ValueError(f"No strategy defined for model class PX")
    def get_BS(self):
        raise ValueError(f"No strategy defined for model class BS")
    def get_BO(self):
        raise ValueError(f"No strategy defined for model class BO")
    def get_LS(self):
        raise ValueError(f"No strategy defined for model class LS")
    def get_LX(self):
        raise ValueError(f"No strategy defined for model class LX")
    def get_LO(self):
        raise ValueError(f"No strategy defined for model class LO")
    def get_LK(self):
        raise ValueError(f"No strategy defined for model class LK")
    def get_TS(self):
        raise ValueError(f"No strategy defined for model class TS")
    def get_TX(self):
        raise ValueError(f"No strategy defined for model class TX")
    def get_TO(self):
        raise ValueError(f"No strategy defined for model class TO")



class SelectInitConds:
    def __init__(self,modnumber, grid_output,eventpath,reference_chi2):
        self.modnumber = modnumber
        self.modelcodes =  ['PS', 'PX', 'BS', 'BO', 'LS', 'LX', 'LO', 'LK', 'TS', 'TX','TO']
        self.grid_output = grid_output
        self.eventpath = eventpath
        self.reference_chi2 = reference_chi2

        run_map = {
        0: self.PS_InitConds,
        1: self.PX_InitConds,
        2: self.BS_InitConds,
        3: self.BO_InitConds,
        4: self.LS_InitConds,
        5: self.LX_InitConds,
        6: self.LO_InitConds,
        7: self.LK_InitConds,
        8: self.TS_InitConds,
        9: self.TX_InitConds,
        10: self.TO_InitConds
        }
        run_map[self.modnumber]()

    def PS_InitConds(self):
        raise ValueError(f"No ICGS strategy defined for model class PS")
    def PX_InitConds(self):
        raise ValueError(f"No ICGS strategy defined for model class PX")
    def BS_InitConds(self):
        raise ValueError(f"No ICGS strategy defined for model class BS")
    def BO_InitConds(self):
        raise ValueError(f"No ICGS strategy defined for model class BO")
    def LS_InitConds(self):
        """
        Finds best initial conditions on grid, with only one model allowed for each q and s.
        Tries to ensure 4 unique mass ratios, but will stop at 20 total initial conditions.
        Mini steps:
        1. Saves 13 best models with unique combination of q and s.
        2a. If < 4 unique mass ratios, check for 3 best models with different mass ratios, then 2, then 2 again.
        In a way there are 4 unique mass ratios.
        2b. If >= 4 unique mass ratios, add best models with unique q,s until 20 initial conditions.
        """



        pspl_thresh = 0
        number_of_best_models_before_q_filter = 13
        number_of_final_models = 20
        number_of_minimum_unique_q = 4
        n_models = self.grid_output.shape[0]
        delta_chi2 = self.reference_chi2 - self.grid_output[:,-1]
        self.grid_output['delta_chi2'] = delta_chi2
        q_values, s_values, alpha_values = self.grid_output.loc[:,'q'], self.grid_output.loc[:,'s'],self.grid_output.loc[:,'alpha']
        # Sort all models by chi2 (ascending)
        argsort_chi2_indices = np.argsort(-delta_chi2)
        delta_chi2_le_0_flag = 0  # flag to set if all models that improve the chi2 have been found.
        best_grid_model_indices = []
        best_grid_model_q_s_tuples = []
        # ----------------------
        # Step 1: select top 13 (number_of_best_models_before_q_filter) unique (q,s)
        # ----------------------
        print('Picking best models. 1 alpha per q,s pair allowed.')
        saved_skipped_indices = []
        last_step1_model_index = 0
        for model_index in range(n_models):
            if len(best_grid_model_indices) == number_of_best_models_before_q_filter:
                last_step1_model_index = model_index - 1  # Just to be true to the variable name
                # print(f'Step 1: {len(best_grid_model_indices)} initial conditions selected')
                break
            current_index = argsort_chi2_indices[model_index]
            q, s = q_values[current_index], s_values[current_index]
            if (q, s) in best_grid_model_q_s_tuples:
                saved_skipped_indices.append(model_index)
                continue
            else:
                if delta_chi2[current_index] > pspl_thresh:
                    best_grid_model_indices.append(current_index)
                    best_grid_model_q_s_tuples.append((q, s))
                else:
                    print("No more models improve the chi2 over the pspl!")
                    delta_chi2_le_0_flag = 1
                    break
        # If we already have all models that improve chi2, print and go to best_model_indices part:
        if delta_chi2_le_0_flag:
            print(f'{len(best_grid_model_indices)} initial conditions selected')

        # ----------------------
        # Step 2: fill remaining slots to reach number_of_final_models (20)
        # ----------------------
        else:
            # Track unique q (mass ratios)
            q_in_init_conds = []
            for qs in best_grid_model_q_s_tuples:
                q_in_init_conds.append(qs[0])
            unique_q_in_init_conds = np.unique(q_in_init_conds)

            # if there are already 4 unique mass ratios, continue adding best (q,s) models
            if len(unique_q_in_init_conds) >= number_of_minimum_unique_q:
                start_i = last_step1_model_index + 1
                for remaining_model_index in range(start_i, n_models):
                    if len(best_grid_model_indices) == number_of_final_models:
                        # print(f'Step 2: {len(best_grid_model_indices)} initial conditions selected')
                        break
                    current_index = argsort_chi2_indices[remaining_model_index]
                    s = s_values[current_index]
                    q = q_values[current_index]
                    if (q, s) in best_grid_model_q_s_tuples:
                        saved_skipped_indices.append(remaining_model_index)
                        continue
                    else:
                        if delta_chi2[current_index] > pspl_thresh:
                            best_grid_model_indices.append(current_index)
                            best_grid_model_q_s_tuples.append((q, s))
                        else:
                            print("No more models improve the chi2 over the pspl!")
                            delta_chi2_le_0_flag = 1
                            break
            # however, if there are less than 4 unique mass ratios, try to add more
            else:
                print("Warning! Less than 4 unique mass ratios in the initial conditions")
                print("Adding more initial conditions")
                # batch_sizes = [3, 2, 2]
                # That's what we have /\ , but writing it more generic:
                number_of_models_that_are_needed = number_of_final_models - number_of_best_models_before_q_filter
                batch_sizes = [number_of_models_that_are_needed // (number_of_minimum_unique_q - 1) +
                               (i < (number_of_models_that_are_needed) % (number_of_minimum_unique_q - 1))
                               for i in range(number_of_minimum_unique_q - 1)]
                start_i = last_step1_model_index + 1
                nstop = number_of_best_models_before_q_filter
                for batch_counter in range(3):
                    nstop += batch_sizes[batch_counter]
                    # print('nstop',nstop)
                    if len(unique_q_in_init_conds) < number_of_minimum_unique_q:
                        for new_model_index in range(start_i, n_models):
                            if len(best_grid_model_indices) == nstop:
                                # If batch is complete, it recalculates again the number of unique mass ratios
                                # and moves the start index to the last model added and breaks to start a new batch
                                q_in_init_conds = []
                                for qs in best_grid_model_q_s_tuples:
                                    q_in_init_conds.append(qs[0])
                                unique_q_in_init_conds = np.unique(q_in_init_conds)
                                start_i = new_model_index
                                break
                            else:
                                # While batch is not complete, keep adding models
                                # if unique q is still not reached, only add models with new mass ratios:
                                current_index = argsort_chi2_indices[new_model_index]
                                q, s = q_values[current_index], s_values[current_index]
                                if ((q, s) in best_grid_model_q_s_tuples) or (q in unique_q_in_init_conds):
                                    saved_skipped_indices.append(new_model_index)
                                    continue
                                else:
                                    if delta_chi2[current_index] > pspl_thresh:
                                        best_grid_model_indices.append(current_index)
                                        best_grid_model_q_s_tuples.append((q, s))
                                        # print(f'Step 3a: {len(best_grid_model_indices)} initial conditions selected')
                                    else:
                                        print("No more models improve the chi2 over the pspl!")
                                        # print(f'Step 3b: {len(best_grid_model_indices)} initial conditions selected')
                                        delta_chi2_le_0_flag = 1
                                        break
                    else:
                        print('Enough unique q added, adding skipped models.')
                        # if minimum unique q is reached, complete with next models from skipped models
                        # calculate how many more models are needed to reach 20
                        number_of_models_needed = number_of_final_models - len(best_grid_model_indices)
                        for skipped_index in range(0, number_of_models_needed):
                            if len(best_grid_model_indices) == nstop:
                                # If batch is complete, it breaks to start a new batch
                                # print(f'Step 3c: {len(best_grid_model_indices)} initial conditions selected')
                                break
                            else:
                                if skipped_index < len(saved_skipped_indices):
                                    saved_index = saved_skipped_indices[skipped_index]
                                else:
                                    saved_index = start_i  # fallback if no more skipped indices

                                current_index = argsort_chi2_indices[saved_index]
                                q, s = q_values[current_index], s_values[current_index]
                                if delta_chi2[current_index] > pspl_thresh:
                                    best_grid_model_indices.append(current_index)
                                    best_grid_model_q_s_tuples.append((q, s))
                                    # print(f'Step 3d: {len(best_grid_model_indices)} initial conditions selected')
                                else:
                                    print("No more models improve the chi2 over the pspl!")
                                    # print(f'Step 3b: {len(best_grid_model_indices)} initial conditions selected')
                                    delta_chi2_le_0_flag = 1
                                    break
        self.best_grid_models = self.grid_output[best_grid_model_indices,:]
    def LX_InitConds(self):
        raise ValueError(f"No ICGS strategy defined for model class LX")
    def LO_InitConds(self):
        raise ValueError(f"No ICGS strategy defined for model class LO")
    def LK_InitConds(self):
        raise ValueError(f"No ICGS strategy defined for model class LK")
    def TS_InitConds(self):
        raise ValueError(f"No ICGS strategy defined for model class TS")
    def TX_InitConds(self):
        raise ValueError(f"No ICGS strategy defined for model class TX")
    def TO_InitConds(self):
        raise ValueError(f"No ICGS strategy defined for model class TO")
