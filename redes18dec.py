#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np


# In[2]:


meta_path = "hmp2_metadata_2018-08-20.csv"


# In[3]:


meta = pd.read_csv(meta_path, low_memory=False)
trans  = pd.read_csv(
    "host_tx_counts.tsv",
    sep="\t",
    index_col=0
)
micro = pd.read_csv(
    "taxonomic_profiles.tsv",
    sep="\t", 
    index_col=False
)


# In[4]:


print("meta:", meta.shape)
print("micro:", micro.shape)
print("trans:", trans.shape)


# In[5]:


meta.head()


# micro.head()

# In[6]:


trans.head()


# In[7]:


micro.head()


# In[8]:


print("Primeras columnas:", micro.columns[:5].tolist())


# In[9]:


trans_T = trans.T

trans_T1= trans.T.reset_index().rename(columns={"index": "External_ID_trans"})


# In[10]:


trans_T1.head()


# In[11]:


# quitar columna OTU
micro = micro.drop(columns=["#OTU ID"])

# poner taxonomy como índice
micro = micro.set_index("taxonomy")
micro.index.name = "External_ID_micro"

micro.shape


# In[12]:


micro_T = micro.T
micro_T.index.name = "External_ID_micro"

micro_T.shape


# In[13]:


micro_T.head()


# In[14]:


# micro_T
d_micro = micro_T.reset_index().duplicated(keep=False).sum()
print("Duplicados exactos en micro_T:", int(d_micro))

# trans_T1
d_trans = trans_T1.reset_index().duplicated(keep=False).sum()
print("Duplicados exactos en trans_T1:", int(d_trans))

# meta
d_meta = meta.duplicated(keep=False).sum()
print("Duplicados exactos en meta:", int(d_meta))


# In[15]:


cols = ['External ID', 'Participant ID', 'week_num', 'visit_num', 'site_sub_coll', 'data_type']

dup_count = meta.duplicated(subset=cols, keep=False).sum()
print("Duplicados en meta por esas columnas:", int(dup_count))


# In[16]:


meta["data_type"].value_counts(dropna=False)


# In[17]:


meta['I tried hard to get to sleep...'].value_counts(dropna=False)


# In[18]:


site_set  = set(meta["External ID"].astype(str).str.strip().dropna())
micro_set = set(micro_T.index.astype(str).str.strip())  # asumiendo index = External_id_micro

matches = site_set & micro_set
matches_list = sorted(matches)

print("Matches:", len(matches))
print("Matches 50 first:", matches_list[:50])
print("site_sub_coll únicos:", len(site_set))
print("eExternal_ID_micro únicos:", len(micro_set))


# In[19]:


dup_count = micro_T.index.duplicated(keep=False).sum()
print("Duplicados en External_id_micro (index):", int(dup_count))


# In[20]:


meta_matches = meta[meta["External ID"].astype(str).str.strip().isin(matches)]
print(meta_matches["data_type"].value_counts(dropna=False))


# In[21]:


meta_set  = set(meta["External ID"].astype(str).str.strip().dropna())
trans_set = set(trans_T1["External_ID_trans"].astype(str).str.strip())

matches_tx = meta_set & trans_set
matches_tx_list = sorted(matches_tx)

print("Matches:", len(matches_tx_list))
print("Matches 50 first:", matches_tx_list[:50])
print("External ID únicos (meta):", len(meta_set))
print("External_ID_trans únicos (trans_T1 index):", len(trans_set))

# data_type de esos matches en meta
meta_matches_tx = meta[meta["External ID"].astype(str).str.strip().isin(matches_tx)]
print("\nConteo data_type en esos matches:")
print(meta_matches_tx["data_type"].value_counts(dropna=False))


# In[22]:


# quedarnos SOLO con host_transcriptomics y biopsy_16S en meta
keep_types = ["host_transcriptomics", "biopsy_16S"]

meta_clean = meta[meta["data_type"].astype(str).str.strip().isin(keep_types)].copy()

print("meta original:", meta.shape)
print("meta_clean:", meta_clean.shape)
print(meta_clean["data_type"].value_counts(dropna=False))


# In[23]:


meta_set  = set(meta_clean["External ID"].astype(str).str.strip().dropna())
trans_set = set(trans_T1["External_ID_trans"].astype(str).str.strip())   

matches_tx = meta_set & trans_set
matches_tx_list = sorted(matches_tx)

print("Matches:", len(matches_tx_list))
print("Matches 50 first:", matches_tx_list[:50])
print("External ID únicos (meta_clean):", len(meta_set))
print("External_ID_trans únicos (trans_T1 index):", len(trans_set))

# data_type de esos matches en meta_clean
meta_matches_tx = meta_clean[meta_clean["External ID"].astype(str).str.strip().isin(matches_tx)]
print("\nConteo data_type en esos matches (meta_clean):")
print(meta_matches_tx["data_type"].astype(str).str.strip().value_counts(dropna=False))


# In[24]:


meta_ids_clean  = set(meta_clean["External ID"].astype(str).str.strip().dropna())
trans_ids       = set(trans_T1["External_ID_trans"].astype(str).str.strip())

missing_in_meta = sorted(trans_ids - meta_ids_clean)

print("IDs en trans_T1 que NO están en meta_clean:", len(missing_in_meta))
print(missing_in_meta)


# In[25]:


meta_ids_all = set(meta["External ID"].astype(str).str.strip().dropna())
missing_in_meta_all = sorted(trans_ids - meta_ids_all)

print("IDs en trans_T1 que NO están ni en meta:", len(missing_in_meta_all))
print(missing_in_meta_all)


# In[26]:


orphans = ['CSMDRVY8', 'CSMDRVY9', 'HSm9JTBX']

trans_T1_clean = trans_T1[~trans_T1["External_ID_trans"].astype(str).str.strip().isin(orphans)].copy()

print("trans_T1 original:", trans_T1.shape)
print("trans_T1 clean (sin huérfanos):", trans_T1_clean.shape)


# In[27]:


meta_clean.to_csv("meta_clean.csv", index=False)
print("Guardado como meta_clean.csv")


# In[30]:


cols = ["Participant ID", "visit_num", "data_type", "week_num"]

dup_count = meta_clean.duplicated(subset=cols, keep=False).sum()
print("Filas duplicadas por Participant ID + visit_num + data_type:", int(dup_count))


dups = meta_clean.loc[meta_clean.duplicated(subset=cols, keep=False), cols].sort_values(cols)
dups.head(50)


# In[33]:


key = ["Participant ID","week_num","visit_num","data_type"]

grp = meta_clean.groupby(key).size().reset_index(name="n")
dups = grp[grp["n"] > 1].sort_values("n", ascending=False)

print("Grupos duplicados (n>1):", dups.shape[0])
print("Top 20 grupos más duplicados:")
print(dups.head(50))


# In[34]:


key = ["Participant ID","week_num","visit_num","data_type"]

grp = meta_clean.groupby(key).size().reset_index(name="n")
print("Grupos duplicados (n>1):", (grp["n"]>1).sum())
print("Filas totales dentro de grupos duplicados:", int(grp.loc[grp["n"]>1, "n"].sum()))


# In[35]:


key = ["Participant ID","week_num","visit_num","data_type"]

diag = (meta_clean.groupby(key)
        .agg(n_rows=("External ID","size"),
             n_external_ids=("External ID","nunique"))
        .reset_index()
        .sort_values(["n_rows","n_external_ids"], ascending=False))

print(diag.head(20))


# In[ ]:




