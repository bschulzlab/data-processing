
# coding: utf-8

# In[1]:


import pandas as pd


# In[2]:


source_file = r"C:\Users\localadmin\Downloads\20190418_17_Schulz_Ruby_IgG_mz600_ETH1_MSnSpectrumInfo.xlsx"


# In[3]:


to_replace = r"C:\Users\localadmin\Downloads\ETH1.xlsx"


# In[4]:


output = r"C:\Users\localadmin\Downloads\ETH1_updated.xlsx"


# In[5]:


source_df = pd.read_excel(source_file, index_col=15)
to_replace_df = pd.read_excel(to_replace)


# In[6]:


for i, r in to_replace_df.iterrows():
    to_replace_df.at[i, "Area"] = source_df.loc[r["First Scan"]]["Area"]


# In[7]:


to_replace_df.to_excel(output, index=False)

