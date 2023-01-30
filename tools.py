import numpy as np
import pandas as pd
import itertools

from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import scipy.stats as st

#import bokeh.io
import bokeh.plotting as plt
#from bokeh.layouts import gridplot
#from bokeh.io import export_png

bokeh.io.output_notebook()

# Importer
def importer(ct_file):
    '''Imports Ct data from Lightcycler machine output of Absolute Quantification:
    Second Derivative Max'''
    
    if ct_file[-3:] == 'csv':
        return pd.read_csv(ct_file, header=1, usecols=['Pos','Name','Cp'],sep='\t')
    
    if ct_file[-3:] == 'txt':
        return pd.read_csv(ct_file, header=1, usecols=['Pos','Name','Cp'],sep='\t')
    
    elif ct_file[-4:] == 'xlsx':
        return pd.read_excel(ct_file, header=1, usecols=['Pos','Name','Cp'],sep='\t')
    
    else:
        raise ValueError('The filetype entered was not recognized. Please enter a .csv, .txt, or .xlsx file.') 
        
# Dilution 
class Dilution:
    def __init__(self,
              with_dil, # samples with dilution curves. Enter as list
              dil_series, # dilution series. Enter as list of integers (i.e. 1:10 dilution = 10)
              dil_rest # the dilution of samples without the dilution series
              ):
        self.with_dil = with_dil
        self.dil_series = dil_series
        self.dil_rest = dil_rest
        
    def update_samples(self,allsamples):
        update_ls = []
        
        for i in allsamples:
            if i in self.with_dil:
                for j in self.dil_series:
                    update_ls.append(i + ' ' + str(j))
            else:
                update_ls.append(i + ' ' + str(self.dil_rest))
                
        return update_ls
    
    
# Sample namer
def namer(ct_file, # file containing Ct values (str)
         primers, # primers, in order (tuple)
         samples, # samples, in order (tuple)
         reps, # technical replicates, [2, 3, or 4]
         config='line', # configuration of technical replicates. choose 'line' or 'square'
         dilution=False # dilution factors for samples if relevant, enter as dilution object 
         ):
    
    '''Imports Lightcycler Abs Quant 2nd Deriv Max data from 384-well plate 
    and labels wells with primers, samples, and dilutions if applicable'''
    
    # Check for valid replicate number
    valid_reps = [2,3,4]
    
    if reps in valid_reps:
        pass
    else:
        raise ValueError('Accepted replicate numbers are 2, 3, and 4.')
    
    
    # Check for invalid conformation
    if reps == 2 and config == 'square':
        ask = input('''You have entered only two technical replicates, 
        but have asserted they are arranged in a square. Would you like to: 
        \n 1) Proceed with line option setting
        \n 2) Cancel analysis and correct replicate number \n''')
        
        if ask == '1':
            config = 'line'
            print('\n Proceeding with line')
        elif ask == '2':
            return print('\n You have chosen to cancel. Please correct your technical replicates.')
        else:
            raise ValueError('Please enter a valid option.')
            
    else:
        pass
    
    # Handle dilutions
    if dilution == False:
        pass
    else:
        samples = dilution.update_samples(samples)
    
    # Read in the data
    ct_data = importer(ct_file)
    
    # Line configuration
    if config == 'line':
        
        # Drop empty wells. If empty wells exist where they shouldn't the program will fail.
        ct_data.dropna(inplace=True)
        
        # Primer list
        primls = []
        
        for p in primers:
            for r in range(reps * len(samples)):
                primls.append(p)
        
        # Sample list
        sampls = []
        
        for s in samples:
            for r in range(reps):
                sampls.append(s)
                
        sampls = sampls * len(primers)
        
        ct_data['Name'] = sampls
        ct_data['Primer'] = primls
     
    # Square configuration
    elif config == 'square':
        # Primers
        primarr = [['nan' for i in range(24)] for j in range(16)]
        
        for i,p in enumerate(primers):
            # First row of primers
            for r in range(len(samples)*2):
                primarr[2*i][r] = p

             # Second row of primers
            if reps == 4:
                for r in range(len(samples)*2):
                    primarr[(2*i)+1][r] = p

            elif reps == 3:
                for si, sv in enumerate(samples):
                    primarr[(2*i)+1][si*2] = p
        
        # Concatenate lists
        primarr = list(itertools.chain.from_iterable(primarr))
        print(len(primarr))
        
        # Make new column and truncate
        ct_data['Primer'] = primarr
        ct_data = ct_data.iloc[:len(primers)*48].copy()
            
        
        # Samples 
        samparr = [['nan' for i in range(24)] for j in range(2)]
        for i,s in enumerate(samples):
            # First row of samples
            for r in range(2):
                samparr[0][i*2] = s
                samparr[0][(i*2)+1]=s

             # Second row of samples
            if reps == 4:
                for r in range(2):
                    samparr[1][i*2] = s
                    samparr[1][(i*2)+1]=s

            elif reps == 3:
                for r in range(2):
                    samparr[1][i*2] = s
        
        samparr = samparr * len(primers)
        
        # Concatenate lists
        samparr = list(itertools.chain.from_iterable(samparr))
        
        # Add to names column
        ct_data['Name'] = samparr
    
    else:
        raise ValueError('Please enter a valid configuration: line or square')
        
    ct_data['NamePrim'] = ct_data['Name'] + ct_data['Primer']
        
    return ct_data[ct_data['Name'] != 'nan']

# Efficiency testing
def efficiency(ct_namer, # output from namer
              dilution, # dilution object
              reps, # number of technical replicates
              usesample=False, # False if using all samples, list of samples if using only some
              returnmodel=False, # Whether or not to output the linear model in full
              ):    
    # Generate dilution column and remove from names
    ct_data = ct_namer.copy()
    name_lsls = [i.split(' ') for i in ct_data['Name']]

    ct_data['Dilution'] = [1/int(i[-1]) for i in name_lsls]
    ct_data['Name'] = [' '.join(i[:-1]) for i in name_lsls]

    # Filter for samples with dilution curves
    if usesample == False:
        ct_data = ct_data[ct_data['Name'].isin(dilution.with_dil)]
    else:
        ct_data = ct_data[ct_data['Name'].isin(usesample)]

    # Calculate log2 dilution
    ct_data.loc[:,'Log2Dil'] = np.log2(ct_data['Dilution'])

    # Initiate dictionaries for efficiency values and plots
    eff_dict = {'Name':[],
               'Primer':[],
               'Efficiency':[],
               'Rsquared':[]}
    
    model_dict = {'Name':[],
                 'Primer':[],
                 'Coefficient':[],
                 'Intercept':[]}

    plot_dict = {}

    # Set linear regression prediction range
    minpred = min(ct_data['Log2Dil'])
    maxpred = max(ct_data['Log2Dil'])

    # Perform linear regression and update plots and efficiencies
    for n in np.unique(ct_data['Name']):
        plot_dict[n] = []
        print(n)
        # Subset by name
        subdf = ct_data[ct_data['Name'] == n]
        # Loop through primers
        for p in np.unique(subdf['Primer']):
            
            primdf = subdf[subdf['Primer'] == p]
            predrange = np.linspace(minpred,maxpred,len(primdf))[::-1]

            # Set up plot for ct values and regression line
            plot = plt.figure(title=' '.join([n,p]), height=400, width=400,
                              x_axis_label = 'Log2Dil', y_axis_label='Ct',)

            plot.circle(primdf['Log2Dil'].values, primdf['Cp'].values)

            # Perform linear regression
            eff_ld = np.array(primdf['Log2Dil']).reshape(-1,1)
            eff_cp = np.array(primdf['Cp'])

            reg = LinearRegression()
            reg.fit(eff_ld, eff_cp)

            bf = reg.predict(predrange.reshape(-1,1))

            # Update plots and calculate efficiency
            plot.line(predrange, bf)
            eff = 2**(-1 / reg.coef_) - 1
            eff = round(eff[0],3)

            rsquared = r2_score(np.array(primdf['Cp']), bf.reshape(len(bf)))

            eff_dict['Name'].append(n)
            eff_dict['Primer'].append(p)
            eff_dict['Efficiency'].append(eff)
            eff_dict['Rsquared'].append(rsquared)
            
            model_dict['Name'].append(n)
            model_dict['Primer'].append(p)
            model_dict['Coefficient'].append(round(reg.coef_[0],3))
            model_dict['Intercept'].append(round(reg.intercept_,3))

            plot_dict[n].append(plot)
            
    # Convert dicts to dfs
    eff_df = pd.DataFrame(eff_dict)
    model_df = pd.DataFrame(model_dict)
    
    # Return efficiency values with or without model
    if returnmodel == True:
        return plot_dict, eff_df, model_df
    else:
        return plot_dict, eff_df
        
# Delta-Delta Ct Analysis
def deltadeltact(ct_data, # output of namer
                 housekeeping, # list of housekeeping primers
                 primers, # list of all primers
                 exp_ctrl, # dictionary of experimental and control samples
                 dilution=False, # dilution factor if relevant, enter as dilution object
                 foldchange=False # whether to output fold change or delta-delta Ct
                ):
    '''Performs delta delta Ct analysis on output of namer program'''
    # Check for dilution
    if dilution == False:
        pass
    else:
        ct_data['Dilutions'] = [int(i.split(' ')[-1]) for i in ct_data['Name']]
        ct_data = ct_data[ct_data['Dilutions'] == dilution.dil_rest]
    
    # Compute averages 
    avgdf = averager(ct_data)
        
    # Check for appropriate housekeeping controls
    if len(housekeeping) == 1:
        ask = input('''Warning: Using only one housekeeping gene severely limits result accuracy.
                    \n Do you want to proceed? [Y/N]''')
        
        if ask == 'Y':
            print('Proceeding with delta delta Ct analysis.')
            
        elif ask == 'N':
            return print('Analysis canceled.')
        
        else:
            raise ValueError('Please enter Y or N to the warning prompt.')
        
    else:
        pass
    
    
    # Batch together housekeeping genes
    names = np.unique(avgdf['Name'])
    
    for n in names:
        subdf = avgdf[avgdf['Name']==n]
        subdf = subdf[subdf['Primer'].isin(housekeeping)]

        errorprop = 0.5 * np.sqrt(np.sum([i**2 for i in subdf['StdCt']]))

        series = pd.Series({'Name':n,
                           'Primer':'housekeeping',
                           'AvgCt':np.mean(subdf['AvgCt']),
                           'StdCt':errorprop
                               })

        avgdf = pd.concat([avgdf,series.to_frame().T],ignore_index=True)
        
    # Calculate delta Ct values (experimental - housekeeping)
    exp_genes = len(primers) - len(housekeeping)
    
    # Make dictionary of empty lists
    dct_dict = {'Name':[],
               'Primer':[],
               'dCt':[],
               'StdErr':[]}
    
    for i, n in enumerate(names):
        # Subset df and make primer loc-able
        subdf = avgdf[avgdf['Name']==n]
        subdf.set_index('Primer',inplace=True)
        
        for p in primers:
            if p in housekeeping:
                pass
            else:
                # Calculate dCt and propagate error
                dct = subdf.loc[p,'AvgCt'] - subdf.loc['housekeeping','AvgCt']
                err = np.sqrt(subdf.loc[p,'StdCt']**2 + subdf.loc['housekeeping','StdCt']**2)
                
                # Update dictionary
                dct_dict['Name'].append(n)
                dct_dict['Primer'].append(p)
                dct_dict['dCt'].append(dct)
                dct_dict['StdErr'].append(err)
                
    # Convert dictionary to dataframe
    dct_df = pd.DataFrame(dct_dict)
    
    # Set index to nameprim for easy location
    dct_df['NamePrim'] = dct_df['Name'] + dct_df['Primer']
    dct_df.set_index('NamePrim',inplace=True)
    
    # Collect the names of the primers used
    dct_prims = np.unique(dct_df['Primer'])
    
    # Make dictionary of empty lists
    ddct_dict = {'Experimental':[],
                'Control':[],
                'Primer':[],
                'ddCt':[],
                'StdErr':[]}
    
    # Calculate delta delta Ct values (experimental dCt - control dCt)
    for e in exp_ctrl:
        # Identify control sample
        c = exp_ctrl[e]
        
        # Loop through primers
        for p in dct_prims:
            # Identify indices
            e_prim = e + p
            c_prim = c + p
            
            # Calculate ddCt and propagate error
            ddct = dct_df.loc[e_prim,'dCt'] - dct_df.loc[c_prim,'dCt']
            err = np.sqrt(dct_df.loc[e_prim,'StdErr']**2 + dct_df.loc[c_prim,'StdErr']**2)
            
            # Update dictionary
            ddct_dict['Experimental'].append(e)
            ddct_dict['Control'].append(c)
            ddct_dict['Primer'].append(p)
            ddct_dict['ddCt'].append(ddct)
            ddct_dict['StdErr'].append(err)
            
    # Convert dictionary to dataframe
    ddct_df = pd.DataFrame(ddct_dict)
    
    # Return the dataframe
    if foldchange == False:
        return ddct_df
    elif foldchange == True:
        print('I have not figured out how to properly propagate error for fold change yet')
        ddct_df['FoldChange'] = 2**(-1 * ddct_df['ddCt'])
        return ddct_df
    else:
        raise ValueError('foldchange should be True or False')