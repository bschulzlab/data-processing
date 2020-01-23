
# coding: utf-8

# The script would read "source_file" as source of "Area" and replace all "Area" in "to_replace" at corresponding "First Scan" values. The output filename and path are defined at "output".

# In[1]:


import pandas as pd


# In[2]:


source_file = r"C:\Users\localadmin\Downloads\20190418_17_Schulz_Ruby_IgG_mz600_ETH1_MSnSpectrumInfo.xlsx"
to_replace = r"C:\Users\localadmin\Downloads\ETH1.xlsx"
output = r"C:\Users\localadmin\Downloads\ETH1_updated.xlsx"


# In[3]:


source_df = pd.read_excel(source_file, index_col=15)
to_replace_df = pd.read_excel(to_replace)


# Iterate through "to_replace_df" and replace "Area" with value at the corresponding "First Scan" from "source_file"

# In[4]:


for i, r in to_replace_df.iterrows():
    try:
        to_replace_df.at[i, "Area"] = source_df.loc[r["First Scan"]]["Area"]
    except:
        print(r["First Scan"] + " does not exist in source file")
    


# In[5]:


to_replace_df.to_excel(output, index=False)

