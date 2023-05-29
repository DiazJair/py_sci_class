# Colección de funciones utiles para procesar los datos de proteomica reportados por Mori et al., 2021.
# Importación de librerías
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Los datos se encuentran en dos datasets distintos:
# De ser necesario, se pueden concatenar con la función concat_datasets().
def concat_datasets():
    """
    Concatena los dos datasets de proteomica.

    Args: 0, la funcion no recibe argumentos. Accedera a los archivos directamente. 
    Returns: Dataframe con ambos datasets concatenados.
    """
    # Obteniendo el directorio actual de trabajo
    cwd = os.getcwd()
    df1 = pd.read_excel(cwd+"/datasets/msb20209536-sup-0009-datasetev8.xlsx", sheet_name="EV8-AbsoluteMassFractions-1")
    df2 = pd.read_excel(cwd+"/datasets/msb20209536-sup-0010-datasetev9.xlsx", sheet_name="EV9-AbsoluteMassFractions-2")
    if df1["Gene name"].equals(df2["Gene name"]):
        # Eliminando las columnas repetidas en el dataset 2.
        df2 = df2.drop(columns=["Gene name", "Gene locus", "Protein ID"])
        mori_df = pd.concat([df1, df2], axis=1)
        print("Los dataframes tienen el mismo orden, se han concatenado.")
        return(mori_df)
    else:
        return("Los DF no tienen el mismo orden, no se pueden concatenar.")
        

def mw_prots(prot_lis, fasta_seq):
    """
    Función que recibe una lista de proteinas en formato tabular, y las secuencias fasta de las proteinas.
    Y devuelve un dataframe con las proteinas de interés y su peso molecular.
    Diseñada para procesar las smart tables creadas por Ecocyc (https://www.ecocyc.org/).

    Args:
        prot_lis: Lista de proteinas en formato tabular.
        fasta_seq: Ruta hacia el archivo en formato fasta.
    Returns:
        df (dataframe): Dataframe con las proteinas de interés.
    """
    # Primera parte: Calcula el peso molecular de las proteinas
    dic_seq = {}
    dic_mw = {}

    with open(fasta_seq, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            dic_seq[record.id] = record.seq
    
    for key, value in dic_seq.items():
        dic_mw[key] = ProteinAnalysis(str(value)).molecular_weight()

    # Segunda parte: Crea un dataframe con las proteinas de interés añadiendo una columna con su peso molecular.
    df = pd.read_csv(prot_lis, sep="\t")
    df["Proteins"] = df["Proteins"].str.replace(" ","_")
    # Conserva solo las columnas tituladas "Proteins" y "Gene"
    df = df[["Proteins", "Gene"]]
    # Añade una columna con el peso molecular de cada proteina
    df["Molecular weight"] = df["Proteins"].map(dic_mw)
    return df

# Función para graficar la aportación de una lista de proteinas al porcentaje total del proteoma.
# La función recibe por parametro una lista de proteinas y un dataframe con los datos de proteomica.
# Devuelve un heatmap con la aportación de cada proteina al porcentaje total del proteoma, en todas las condiciones presentes en el dataframe.
# Adicionalmente, se le puede pasar por parametro una lista de condiciones para filtrar el dataframe.
def prot_heatmap(prot_lis, df, conditions=None):
    """
    Función para graficar la aportación de una lista de proteinas al porcentaje total del proteoma.
    La función recibe por parametro una lista de proteinas y un dataframe con los datos de proteomica.
    Devuelve un heatmap con la aportación de cada proteina al porcentaje total del proteoma, en todas las condiciones presentes en el dataframe.
    Adicionalmente, se le puede pasar por parametro una lista de condiciones para filtrar el dataframe.
    Util para obtener un preview rápido de la aportación de las proteinas de interés al proteoma.

    Args:
        prot_lis: Lista de proteinas de interés.
        df: Dataframe con los datos de proteomica.
        conditions: Lista de condiciones para filtrar el dataframe.
    Returns:
        heatmap: Heatmap con la aportación de cada proteina al porcentaje total del proteoma.
    """
    df_list = []
    # Elimina las columnas que no son de interes.
    df = df.drop(columns=["Gene locus", "Protein ID"])
    for prot in prot_lis:
        df_list.append(df[df["Gene name"] == prot])
    # Concatena los dataframes de cada proteina en uno solo.
    df_hm = pd.concat(df_list)
    df_hm = df_hm.set_index("Gene name")
    # Filtra el dataframe por las condiciones de interes, si es que estas existen.
    if conditions == None:
        heatmap = sns.heatmap(df_hm*100, annot=False, cmap="Blues", linewidths=.5)
    else:
        df_hm = df_hm.filter(conditions)
        heatmap = sns.heatmap(df_hm*100, annot=False, cmap="Blues", linewidths=.5)
    return heatmap