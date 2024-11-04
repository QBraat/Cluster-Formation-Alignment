
from cc3d.core.PySteppables import *
from pathlib import Path
import random
import os
import glob
from time import localtime
import statistics

Np = int({{XDIM}}/10)
kappa = {{KAPPA}}
tau = {{TAU}}
frac = {{FRACTIONACTIVE}}
gMM = {{GMM}}



class ClusterCountingAlgorithmSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)
        

    def start(self):
        global write_clusters, write_cells
        
        global write_single_export, mean_cluster_size_time, order_parameter_time, write_St, write_Smax
       
        mean_cluster_size_time = []
        order_parameter_time = []
        

        pixels = int(10)
        for i in range(Np):
            for j in range(Np):
                    self.cell_field[i*pixels:(i+1)*pixels, j*pixels:(j+1)*pixels,0] = self.new_cell(self.EPITHELIAL)
         
        Number_Mesenchymal = int(frac*Np*Np)
        
        while Number_Mesenchymal > 0: 
            for cell in self.cell_list: 
                if self.dim.x/5 < cell.xCOM < 4*self.dim.x/5 and self.dim.y/5 < cell.yCOM < 4*self.dim.y/5 and cell.type  != self.MESENCHYMAL:
                    cell.type = self.MESENCHYMAL
                    Number_Mesenchymal -= 1
                    break 
                    
        volume_per_cell = int(self.dim.x*self.dim.y/(Np*Np))
            
        for cell in self.cell_list_by_type(self.EPITHELIAL,self.MESENCHYMAL):
            
            if cell.type == self.EPITHELIAL:
                cell.targetVolume = volume_per_cell
                cell.lambdaVolume = 2.0
            else:
                cell.targetVolume = volume_per_cell
                cell.lambdaVolume = 2.0
               
                
                
                
        write_single_export, file_path = self.open_file_in_simulation_output_folder('single-export_N%03d_kappa%04d_tau%04d_frac%03d_10000gMM%05d.dat'%(int(Np),int(kappa),int(tau),int(frac*100),int(gMM*10000)), mode='w')
        if write_single_export is None:
            return   
        print('steadyS steadyS_std steadyP steadyP_std',file=write_single_export)   
           
           
        write_St, file_path = self.open_file_in_simulation_output_folder('St-single-export_N%03d_kappa%04d_tau%04d_frac%03d_10000gamma%05d.dat'%(int(Np),int(kappa),int(tau),int(frac*100),int(gMM*10000)), mode='w')
        if write_St is None:
            return   
                       
        print('mcs St St_std', file=write_St)


        write_Smax, file_path = self.open_file_in_simulation_output_folder('S-max-single-export_N%03d_kappa%04d_tau%04d_frac%03d_10000gamma%05d.dat'%(int(Np),int(kappa),int(tau),int(frac*100),int(gMM*10000)), mode='w')
        if write_Smax is None:
            return   
                       
        print('mcs Smax', file=write_Smax)
        
                
                
        write_metadata, file_path = self.open_file_in_simulation_output_folder('metadata_N%03d_kappa%04d_tau%04d_frac%03d_10000gMM%05d.dat'%(int(Np),int(kappa),int(tau),int(frac*100),int(gMM*10000)), mode='w')
        if write_metadata is None:
            return   
            
        print('Manually set parameters  ',file=write_metadata)
        print('kappa:                   ', kappa,file=write_metadata)
        print('tau:                     ', tau,file=write_metadata)
        print('number_of_particles:     ', Np**2,file=write_metadata)
        print('frac of active cells:    ', frac,file=write_metadata)
        print('gEE:                     ', 0.0,file=write_metadata)
        print('gEM:                     ', 0.0,file=write_metadata)
        print('gMM:                     ', gMM,file=write_metadata)
        
        print('\nCell types',file=write_metadata)
        print('Medium:                  type = 0',file=write_metadata)
        print('Epithelial:              type = 1',file=write_metadata)
        print('Mesenchymal:             type = 2',file=write_metadata)
        
        temp = self.get_xml_element('temp')
        JMedE = self.get_xml_element('JMedE')
        JMedM = self.get_xml_element('JMedM')
        JEE = self.get_xml_element('JEE')
        JEM = self.get_xml_element('JEM')
        JMM = self.get_xml_element('JMM')
        neighorder= self.get_xml_element('neighorder')
        max_mcs = self.get_xml_element('numbersteps')
        
        
        print('\nXML parameters',file=write_metadata)
        print('box-size x:              ', self.dim.x,file=write_metadata)
        print('box-size y:              ', self.dim.y,file=write_metadata)
        print('box-size z:              ', self.dim.z,file=write_metadata)
        print('J_MedE:                  ',float(JMedE.cdata),file=write_metadata) 
        print('J_MedM:                  ', float(JMedM.cdata),file=write_metadata) 
        print('J_EE:                    ', float(JEE.cdata),file=write_metadata) 
        print('J_EM:                    ', float(JEM.cdata),file=write_metadata) 
        print('J_MM:                    ', float(JMM.cdata),file=write_metadata) 
        print('neighbor order:          ', float(neighorder.cdata),file=write_metadata) 
        print('total mcs:               ', float(max_mcs.cdata),file=write_metadata) 

        print('temperature:             ', float(temp.cdata),'\n',file=write_metadata) 

                
        print('time stamp:              {}h{}m{}s'.format(localtime()[3],localtime()[4],localtime()[5]),file=write_metadata)        
        print('day stamp:               {}d{}m{}y'.format(localtime()[2],localtime()[1],localtime()[0]),file=write_metadata)   

        

    def step(self,mcs):
        if mcs == 2500:
            for cell in self.cell_list:
                cell.dict['angle_pol_now'] = np.pi*random.uniform(-1,1)
                
                cell.dict['angle_pol_update'] = cell.dict['angle_pol_now']

        if mcs > 2500:
            for cell in self.cell_list:
                Pix = np.cos(cell.dict['angle_pol_now'])
                Piy = np.sin(cell.dict['angle_pol_now'])
                
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor:
                        if neighbor.type == self.MESENCHYMAL and cell.type == self.MESENCHYMAL:
                            Pix += gMM*np.cos(neighbor.dict['angle_pol_now'])
                            Piy += gMM*np.sin(neighbor.dict['angle_pol_now'])
                
                theta = np.arctan2(Piy,Pix) + 1 * (np.sqrt(2./(tau))*np.random.normal(loc = 0.0, scale=1.0))
                cell.dict['angle_pol_update'] = theta
        
            for cell in self.cell_list:
                if cell.type == self.MESENCHYMAL:
                    cell.lambdaVecX = -kappa*np.cos(cell.dict['angle_pol_update'])
                    cell.lambdaVecY = -kappa*np.sin(cell.dict['angle_pol_update'])
                    cell.lambdaVecZ = 0
                    cell.dict['angle_pol_now'] = cell.dict['angle_pol_update']
                if cell.type == self.EPITHELIAL:
                    cell.dict['angle_pol_now'] = cell.dict['angle_pol_update']




        if mcs>=2500 and mcs%200 == 0:     

            ## creating neighborlists
            cell_and_its_neighbors = []
            mesenchymal_cell_id = []
            for cell in self.cell_list_by_type(self.MESENCHYMAL):
                mesenchymal_cell_id.append(cell.id)
            
                neighbor_list = []
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor and neighbor.type == self.MESENCHYMAL:
                        neighbor_list.append(neighbor.id)
                        
                        
                cell_and_its_neighbors.append([cell.id,neighbor_list])
                
                
            cluster_check_list = []
            cluster_id = 0
            cluster_size = 0
            
            cluster_distribution = []
            cluster_id_per_cell = []
            cell_id_already_checked = []
            
            while len(mesenchymal_cell_id)>0:
                
                ### check if there is any element in the cluster_check list. If not, start new cluster
                if len(cluster_check_list) == 0:
                    
                    cluster_id +=1
                    cellid = mesenchymal_cell_id[0] ## choose a new cell to check the cluster size of the cluster it is part of
                    cluster_size = 1
                    
                    ### check neighbors and fill cluster_check_list if applicable
                    for i in range(len(cell_and_its_neighbors)):
                        if cell_and_its_neighbors[i][0] == cellid:
                            cluster_check_list += cell_and_its_neighbors[i][1] ## add neighbors, i.e. cells in the cluster, to the cluster_list that needs to be checked
                           
                    ## keep track of cells that have been considered already
                    mesenchymal_cell_id.remove(cellid)
                    cell_id_already_checked.append(cellid)       
                    
                    if len(cluster_check_list) == 0:
                        ## if a cluster has no more new neighbors, all cells in the cluster have been included and the cluster is closed
                        cluster_distribution.append([cluster_id, cluster_size])                
                    
                    cluster_id_per_cell.append([cellid,cluster_id])
                    
                else: 
                    cellid = cluster_check_list[0] ## choose a cell that is part of an existing cluster
                    cluster_size += 1
                    
                    ### check neighbors and fill cluster_check_list if applicable
                    for i in range(len(cell_and_its_neighbors)):
                        if cell_and_its_neighbors[i][0] == cellid:
                            cluster_check_list += cell_and_its_neighbors[i][1]

                    ## keep track of cells that have been considered already
                    mesenchymal_cell_id.remove(cellid)
                    cell_id_already_checked.append(cellid)                

                    cluster_check_list = list(set(cluster_check_list))
                    for checked_id in cell_id_already_checked:
                        ## remove all cells that have been checked 
                        if checked_id in cluster_check_list:
                            cluster_check_list.remove(checked_id)
                    
                    if len(cluster_check_list) == 0:
                        ## if a cluster has no more new neighbors, all cells in the cluster have been included and the cluster is closed
                        cluster_distribution.append([cluster_id, cluster_size])  
                    cluster_id_per_cell.append([cellid,cluster_id])

        if mcs>=2500 and mcs%200 == 0:     
            mean_cluster_size_time.append(np.mean(np.array(cluster_distribution)[:,1]))
            px = 0
            py = 0
            
            for cell in self.cell_list_by_type(self.MESENCHYMAL): 
                px += np.cos(cell.dict['angle_pol_now'])
                py += np.sin(cell.dict['angle_pol_now'])
                
            Porder = 1/(frac*Np*Np)*np.sqrt(px**2+py**2)
                
            order_parameter_time.append(Porder)    
     
            print('{} {} {}'.format(mcs,np.mean(np.array(cluster_distribution)[:,1]), np.std(np.array(cluster_distribution)[:,1])) ,file=write_St)
            
            print('{} {}'.format(mcs,max(np.array(cluster_distribution)[:,1])) ,file=write_Smax)

            if mcs == 60000:
                idx = len(mean_cluster_size_time)
            
            
            
    def finish(self):
        """
        Finish Function is called after the last MCS
        """
        mean_cluster_steady = statistics.mean(np.array(mean_cluster_size_time[idx-1:]))
        std_cluster_steady = statistics.stdev(np.array(mean_cluster_size_time[idx-1:]))
        order_parameter_steady = statistics.mean(np.array(order_parameter_time[idx-1:]))
        std_order_parameter_steady = statistics.stdev(np.array(order_parameter_time[idx-1:]))
        
        print('{} {} {} {}'.format(mean_cluster_steady, std_cluster_steady,order_parameter_steady,std_order_parameter_steady),file=write_single_export)   

        write_single_export.close()
        write_St.close()
        write_Smax.close()

    def on_stop(self):
        # this gets called each time user stops simulation
        return