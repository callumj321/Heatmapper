# Import
import MDAnalysis as mda
from MDAnalysis import analysis
from MDAnalysis.analysis import contacts
import multiprocessing
from multiprocessing import Pool
from multiprocessing import cpu_count
from functools import partial
import numpy as np
import pandas as pd
import seaborn as sns
import nglview as nv
import os

# quotes = ['"Lets-a go" - Chris Pratt','"I look very much forward to showing my financials, because they are huge." - Donald Trump','"I did not copy Gromacs" - Some rad guy, probably',
# '"ff at 15" - Abraham Lincoln']

# Define class
class traj_analysis:
    def __init__(self,pdb_path,xtc_path,probe_sel='protein',analyte_sel='protein'):
        '''
        Initialiase the class.
        pdb_path (str)  : String identifying the file to use for the pdb by a local path.
        xtc_path (str)  : String identifying the file to use for the xtc by a local path.
        probe (str)     : The selection for the probe, using MDAnalysis selection language.
        analyte (str)   : The selection for the analyte, using MDAnalysis selection language.
        '''
        self.pdb_path = pdb_path
        self.xtc_path = xtc_path
        # Not sure if the filename without the path will be helpful, but maybe
        self.pdb_filename = os.path.split(pdb_path[-1])
        self.xtc_filename = os.path.split(xtc_path[-1])
        # By default the probe and analyte are 'protein'
        # This gives protein-protein interactions
        self.probe_sel = probe_sel
        self.analyte_sel = analyte_sel
        self.n_jobs = cpu_count() # <--- CPUs available for multithreading

    
    def def_load_in_uni(self):
        '''
        A function to load in an MDA universe if one isn't already loaded.
        '''
        try:
            print(self.mda_universe)
            print('...universe already loaded...')
        except:
            self.mda_universe = mda.Universe(self.pdb_path,self.xtc_path)
    
    def frame_to_time(self,frame_in):
        '''
        A function convert a frame to time (ps) of simulation for providing input to contact analysis.
        '''
        self.def_load_in_uni()
        self.time = frame_in*self.mda_universe.trajectory.dt
        return self.time
    
    def time_to_frame(self,time_in):
        '''
        A function convert a time (ps) to frame of simulation for providing input to contact analysis.
        '''
        self.def_load_in_uni()
        self.frame = time_in/self.mda_universe.trajectory.dt
        return self.frame
    
    def pre_cont_pro(self,start,stop,skip):
        '''
        A function to report useful information on proposed imput parameters for a system before it is analysed.
        This saves the used from having to experiment with run time so much.
        start (int) : The frame in simulation to start analysis on.
        stop (int)  : The frame in simulation to stop analysis on.
        skip (int)  : The frames to skip between calculations. This can save run time.
        '''
        print('Frames: '+str(int(len(self.mda_universe.trajectory[start:stop:skip]))))
        print('Cores: '+str(int(self.n_jobs)))
        print('Passes: '+str(int(len(self.mda_universe.trajectory[start:stop:skip])/self.n_jobs)))
        print('Chains: '+str(int(len(self.analyte_segids))))
        print('Residues per chain: '+str(int(self.analyte_resids)))
    
    def get_analyte_info(self):
        '''
        A function to provide information about the analyte such as protein residue, atoms, etc.
        '''
        self.def_load_in_uni()
        try:
            print(self.analyte_loaded)
            print('...analyte already selected...')
        except:
            self.analyte_loaded = self.mda_universe.select_atoms(self.analyte_sel)
        self.analyte_segids = np.unique(self.analyte_loaded.residues.segids)
        self.analyte_resids = np.unique(self.analyte_loaded.residues.resids)

    def cont_per_frame(self,frame_index,cont_dist,carbon,segid='A'): # <--- The actual function which is executed on each CPU
        '''
        A function to perform contact analysis per frame to allow for multithreading.
        This function assumes to first chain only but the actual value is given from the contact function.
        frame_index (int)   : The frame for which the node is calculating.
        cont_dist (float)   : The distance (Angstroms) over which to consider a contact.
        carbon (boolean)    : Whether to include carbons or not. This decreases accuracy but can save on time.
        segid (str)         : The segid of the protein to calculate for.
        '''
        # print('Reached frame '+str(frame_index)+'           ',end='\r')
        print(str(frame_index),end='\r')
        self.mda_universe.trajectory[frame_index] # <--- Select the frame for analysis
        residue_contacts=[] # <--- Create an empty array to store the contacts per residue in
        for resid_iter in self.analyte_resids: # <--- Iterate throught resids
            if carbon == True:
                group_A = self.mda_universe.select_atoms('('+str(self.analyte_sel)+') and (resid '+str(resid_iter)+')') # <--- Select the correct resid
                group_B = self.mda_universe.select_atoms('(around '+str(2*cont_dist)+' resid '+str(resid_iter)+') and ('+str(self.probe_sel)+')') # <--- Select the analyte around the resid
            elif carbon == False:
                group_A = self.mda_universe.select_atoms('('+str(self.analyte_sel)+') and (resid '+str(resid_iter)+') and (not name C)') # <--- Select the correct resid
                group_B = self.mda_universe.select_atoms('(around '+str(2*cont_dist)+' resid '+str(resid_iter)+') and ('+str(self.probe_sel)+') and (not name C)') # <--- Select the analyte around the resid
            distances = mda.analysis.distances.distance_array(group_A.positions, group_B.positions) # <--- Get distances
            contact_count = np.count_nonzero(distances <= cont_dist) # <--- Count the number of distances under the cutoff
            residue_contacts.append(contact_count) # <--- Add the number of contacts for that residus
        print('b')
        df_out = pd.read_csv('tmp\\contact_analysis_data_'+str(segid)+'.csv')
        df_out[str(frame_index)] = residue_contacts
        df_out.to_csv('tmp\\contact_analysis_data_'+str(segid)+'.csv',index=False)
        
    def cont_pro(self,cont_dist=3.3,carbon=True,start=0,stop=-1,skip=1):
        '''
        A function to perform the main contact analysis. The "big papa" of the contained functions, if you will...
        cont_dist (float)   : The distance (Angstroms) over which to consider a contact.
        carbon (boolean)    : Whether to include carbons or not. This decreases accuracy but can save on time.
        start (int)         : The frame in simulation to start analysis on.
        stop (int)          : The frame in simulation to stop analysis on.
        skip (int)          : The frames to skip between calculations. This can save run time.
        '''
        self.def_load_in_uni()
        self.get_analyte_info()
        print('Segids to analyse: ',str(int(len(self.analyte_segids))))
        print('Resids to analyse: ',str(int(len(self.analyte_resids))))
        segid_iter = 0
        for segid in self.analyte_segids: # <--- Iterate through chains
            df_out = pd.DataFrame({'Resid':self.analyte_resids})
            df_out.to_csv('tmp\\contact_analysis_data_'+str(segid)+'.csv',index=False)
            print('Chain: '+str(segid),end='\n') # <--- Update console
            run_per_frame = partial(self.cont_per_frame, cont_dist=cont_dist,carbon=carbon,segid=segid) # <--- Set the per frame function
            print(self.mda_universe.trajectory.n_frames)
            frame_values = np.arange(self.mda_universe.trajectory.n_frames) # <--- Get all frames
            self.analysis_frame_values = frame_values[start:stop:skip] # <--- Select every nth frame for reducing run time
            print(str(len(self.analysis_frame_values))+' frames to analyse on '+str(self.n_jobs)+' cores.',end='\n') # <--- Update console
            print('First frame: '+str(self.analysis_frame_values[0])+' - Final frame: '+str(self.analysis_frame_values[-1]))
            print('Working on it...')
            with Pool(self.n_jobs) as worker_pool: # <--- Create a pool of CPUs to use
                worker_pool.map(run_per_frame, self.analysis_frame_values) # <--- Run the per frame function on a select CPU
            print('here')
            segid_iter+=1
            print('\n')

    def gen_map(self,df_data,res_op=1,res_rep='spacefill',res_rad=0,BB_only=True):
        '''
        A function to generate a heatmap from the data, currently specifically for a protein.
        df_data (dataframe) : A pandas dataframe of the contact data, as output by the cont_pro function.
        res_op (float, 0-1) : The opacity of residues in the output heatmap.
        res_rep (str)       : The NGLview representation to apply to residues.
        res_rad (int)       : The radius of residues to show.
        BB_only (boolean)   : Whether to show the backbone of the protein or not.
        '''
        # Standardise df so that the resids are the index.
        print('Standardising input')
        if 'resids' in df_data.columns.values:
            df_data = df_data.set_index('resids')
        print('Getting contact totals')
        conts_per_resid = []
        for resid in self.analyte_resids:
            conts_per_resid.append(df_data.query('index == @resid').values.sum())
        # return conts_per_resid
        print('Got contact totals',end='\r')
        self.min_conts = min(conts_per_resid)
        self.max_conts = max(conts_per_resid)
        conts_per_resid_minus_min = conts_per_resid-self.min_conts # <--- Scale starts at 0
        conts_per_resid_std = ((conts_per_resid_minus_min/max(conts_per_resid_minus_min))*100).astype(int) # <--- 100 colour values
        sequential_colors = sns.color_palette("flare", 101).as_hex()
        act = self.mda_universe.select_atoms(self.analyte_sel)
        self.heatmap = nv.show_mdanalysis(act,default_representation=True)
        self.heatmap.add_representation('trace', selection='all', color='grey', opacity=1)
        resid_iter = 0
        if res_rad == 0:
            if BB_only == True:
                for resid in self.analyte_resids:
                    self.heatmap.add_representation(res_rep,selection='backbone and '+str(resid),color=sequential_colors[conts_per_resid_std[resid_iter]],opacity=res_op)
                    resid_iter+=1
            else:
                for resid in self.analyte_resids:
                    self.heatmap.add_representation(res_rep,selection=str(resid),color=sequential_colors[conts_per_resid_std[resid_iter]],opacity=res_op)
                    resid_iter+=1
        else:
            if BB_only == True:
                for resid in self.analyte_resids:
                    self.heatmap.add_representation(res_rep,selection='backbone and '+str(resid),color=sequential_colors[conts_per_resid_std[resid_iter]],opacity=res_op,radius=res_rad)
                    resid_iter+=1
            else:
                for resid in self.analyte_resids:
                    self.heatmap.add_representation(res_rep,selection=str(resid),color=sequential_colors[conts_per_resid_std[resid_iter]],opacity=res_op,radius=res_rad)
                    resid_iter+=1
        self.heatmap._view_height = '400px'
        print('Done. Run \nobj.heatmap \nin jupyter to view. Run \nnv.write_html("index.html",[obj.heatmap]) \nTo save as an html you can open later.')
        return self.heatmap