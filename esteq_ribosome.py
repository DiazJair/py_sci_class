## Importando las librerias necesarias
import pandas as pd
import os
import random
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import prot_functs as pf

## Imprimiendo un mensaje de bienvenida
print("Bienvenido al programa de Jair Ernesto Díaz Ramírez.")
print("Iniciando el programa...")

## Definiendo la ruta de trabajo
cwd = os.getcwd() # Guardando el directorio actual en una variable.

## Creando una carpeta en la que se guardaran las imagenes generadas.

# Ruta de la carpeta
carpeta = cwd + "/imagenes"

# Verificar si la carpeta ya existe
if not os.path.exists(carpeta):
    # Crear la carpeta si no existe
    os.mkdir(carpeta)
    print("Carpeta para imagenes creada")
else:
    print("La carpeta 'imagenes' ya existe")

# Obeteniendo los datos de los genes asociados a los ribosomas a partir de los archivos de entrada
ribo_prots = pf.mw_prots(prot_lis=cwd+"/inputs/smart_table.txt", fasta_seq=cwd+"/inputs/proteinas_interes.fasta")

# Obteniendo los datos de proteomica
mori_set = pf.concat_datasets()

# Buscando aquellos genes ribosomales que se encuentren en los datos de proteomica.
ribo_mori_set = mori_set[mori_set["Gene name"].isin(ribo_prots["Gene"])]

# Añadiendo una columna con la subunidad a la que corresponde cada proteina
ribo_mori_set = ribo_mori_set.copy()
ribo_mori_set["Subunit"] = ""

for gene in ribo_mori_set["Gene name"]:
    if ribo_prots[ribo_prots["Gene"] == gene]["Proteins"].str.contains("30S").any():
        ribo_mori_set.loc[ribo_mori_set["Gene name"] == gene, "Subunit"] = "30S"
    elif ribo_prots[ribo_prots["Gene"] == gene]["Proteins"].str.contains("50S").any():
        ribo_mori_set.loc[ribo_mori_set["Gene name"] == gene, "Subunit"] = "50S"
    else:
        ribo_mori_set.loc[ribo_mori_set["Gene name"] == gene, "Subunit"] = "other"

# Añadiendo una columna al dataframe que contenga el peso molecular de cada proteína.
ribo_mori_set["MW"] = ribo_mori_set["Gene name"].map(ribo_prots.set_index("Gene")["Molecular weight"])

# Asignando el nombre del gen como índice del dataframe.
ribo_mori_set.set_index("Gene name", inplace=True)

# Eliminando las filas que correspondan a subunit = "other"
ribo_mori_set = ribo_mori_set[ribo_mori_set["Subunit"] != "other"]

# Eliminando las columnas que no se emplearan
ribo_mori_set.drop(columns=["Gene locus", "Protein ID"], axis=1, inplace=True)

# Calculando la concentración de cada proteína
# Es necesario convertir el molecular weight dado en Daltones, a fentogramos. Para ello se multiplica el valor por 1.6605402e-9.

ribo_mori_set["MW"] = ribo_mori_set["MW"].apply(lambda x: x * 1.6605402e-9)
ribo_mori_set_norm = ribo_mori_set.drop(["MW", "Subunit"], axis=1).apply(lambda x: (x*135)/ribo_mori_set["MW"])

# Añadiendo la columna Subunit al nuevo dataframe.
ribo_mori_set_norm["Subunit"] = ribo_mori_set["Subunit"]

# Graficando el dataframe en un heatmap.
fig, ax = plt.subplots(figsize=(15, 15))
sns.heatmap(data = ribo_mori_set_norm.drop(["Subunit"], axis=1), cmap='Blues', annot=False, linewidths=.5, ax=ax)
ax.set_xlabel("Condicion de crecimiento", fontsize=15)
ax.set_ylabel("Proteina ribosomal", fontsize=15)
plt.suptitle("Abundancia de proteínas ribosomales en E. coli \nNúmero de copias por micrometro cubico", fontsize=20)
plt.savefig(cwd+"/imagenes"+"/ribo_mori_set_norm.png", dpi=300, bbox_inches='tight')

# Eliminando fila correspondiente a ykgO. Se comprobó que en todos los casos, la abundancia de esta proteina es 0.
ribo_mori_set_norm = ribo_mori_set_norm[ribo_mori_set_norm.index != "ykgO"]

# Observando la dispersión de los datos de 5 condiciones al azar en un boxplot.
# Seleccionando 15 condiciones al azar.
condiciones = ribo_mori_set_norm.columns[:-1]
condiciones_5 = random.sample(list(condiciones),15)
# Creando un dataframe con las condiciones seleccionadas.
ribo_mori_set_norm_5 = ribo_mori_set_norm[condiciones_5]
# Graficando el dataframe en un boxplot.
fig, ax = plt.subplots(figsize=(10, 5))
sns.boxplot(data=ribo_mori_set_norm_5, ax=ax)
ax.set_xlabel("Condicion de crecimiento", fontsize=12)
ax.set_ylabel("Número de copias por micrometro cubico", fontsize=12)
plt.suptitle("Abundancia de proteínas ribosomales en E. coli \n15 Condiciones seleccionadas al azar", fontsize=14)
plt.savefig(cwd+"/imagenes"+"/dispersion_15_cond.png", dpi=300, bbox_inches='tight')

# Evaluando la diferencia entre las subunidades 30S y 50S.
# Seleccionando las proteinas 30S y 50S del dataframe ribo_mori_set_norm.
ribo_mori_set_norm_30S = ribo_mori_set_norm[ribo_mori_set_norm["Subunit"] == "30S"].drop(["Subunit"], axis=1)
ribo_mori_set_norm_50S = ribo_mori_set_norm[ribo_mori_set_norm["Subunit"] == "50S"].drop(["Subunit"], axis=1)
# Aplicando el test de Welch.
sp.stats.ttest_ind(ribo_mori_set_norm_30S, ribo_mori_set_norm_50S, equal_var=False)
# Graficando los resultados del test de Welch.
fig, ax = plt.subplots(figsize=(10, 5))
sns.boxplot(data=[ribo_mori_set_norm_30S, ribo_mori_set_norm_50S], ax=ax)
ax.set_xlabel("Subunidad ribosomal", fontsize=12)
ax.set_ylabel("Número de copias por micrometro cubico", fontsize=12)
plt.suptitle("Diferencia en la abundancia de proteínas ribosomales por subunidad", fontsize=14)
ax.set_xticklabels(["30S", "50S"])
plt.savefig(cwd+"/imagenes"+"/dispersion_30S_50S.png", dpi=300, bbox_inches='tight')

# Aplicando el test de Welch.
welch_results = sp.stats.ttest_ind(ribo_mori_set_norm_30S, ribo_mori_set_norm_50S, equal_var=False)
# accediendo a los valores de p-value de welch_results y guardandolos en una variable.
p_values = welch_results[1]

# Nivel de significancia (alfa)
alfa = 0.05

# Comparar cada valor p con alfa
medias_iguales = 0
medias_diferentes = 0
for p in p_values:
    if p < alfa:
        medias_diferentes+=1
        #print("Se rechaza la hipótesis nula. Las medias son diferentes.")
        
    else:
        medias_iguales+=1
        #print("No se puede rechazar la hipótesis nula. Las medias son iguales.")
        
#print("Los valores de p que son menores que alfa son: " + str(medias_diferentes))
#print("Los valores de p que son mayores que alfa son: " + str(medias_iguales))

# Agrupa los genes por subunidad y calcula la media de cada grupo.
r_media_n_set = ribo_mori_set_norm.groupby("Subunit").mean()
# Creando una nueva fila que contenga la relación entre las subunidades 30S y 50S.
r_media_n_set.loc["30S/50S"] = r_media_n_set.loc["30S"] / r_media_n_set.loc["50S"]

# Graficando la fila 30S/50S del dataframe r_media en un barplot.
fig, ax = plt.subplots(figsize=(5, 12))
sns.barplot(data=r_media_n_set.loc["30S/50S"].to_frame().T, orient="horizontal", color="lightblue", ax=ax)
# Trazando una linea vertical en el valor 1, que es lo que se espera.
plt.axvline(x=1, color="red", linestyle="--")
plt.title("30S/50S ratio - Mori et al. (2019)")
plt.xlabel("30S/50S ratio")
plt.ylabel("Condición de crecimiento")
plt.savefig(cwd+"/imagenes"+"/30S_50S_ratio.png", dpi=300, bbox_inches='tight')


ratio_mori = r_media_n_set.loc["30S/50S"]
# Calculando la desviación estándar
std_mori = ratio_mori.std()
# Calculando la varianza de ratio_mori.
var_mori = ratio_mori.var()
# Calculando el coeficiente de variación de ratio_mori.
cv_mori = std_mori / ratio_mori.mean()
# Calculando los percentiles 25, 50 y 75 de ratio_mori.
q1_mori = ratio_mori.quantile(0.25)
q2_mori = ratio_mori.quantile(0.5)
q3_mori = ratio_mori.quantile(0.75)

# Obteniendo los outliers de ratio_mori.
ratio_mori_posible_outliers = ratio_mori[(ratio_mori < ratio_mori.mean() - 1 * ratio_mori.std()) | (ratio_mori > ratio_mori.mean() + 1 * ratio_mori.std())]
list_pos_outliers = ratio_mori_posible_outliers.index.to_list()
set_pos_outliers = ribo_mori_set_norm.loc[:,list_pos_outliers]
set_pos_outliers = set_pos_outliers.reset_index()

outliers = pd.DataFrame(columns=["Gene name"])
outliers["Gene name"] = set_pos_outliers["Gene name"]

for col in set_pos_outliers.columns[1:]:
    mean = set_pos_outliers[col].mean() # Calculando la media de la condición.
    std = set_pos_outliers[col].std() # Calculando la desviación estándar de la condición.
    outliers[col] = set_pos_outliers[col].apply(lambda x: True if x > mean + 2 * std or x < mean - 2 * std else False) # Evaluando cada valor de la condición.

# Por cada fila del dataframe outliers, se obtiene el nombre de la proteina dada por "Gene name" luego se obtiene el conteo de True.
# El conteo de True es el número de condiciones en las que la proteina esta por encima o por debajo de la media más o menos 2 desviaciones estándar.
# Crea un diccionario con el nombre de la proteina y el conteo de True.
outliers_dict = {}
for index, row in outliers.iterrows():
    outliers_dict[row["Gene name"]] = row.iloc[1:].value_counts().get(True, 0)

# Eliminando aquellas claves del diccionario outliers_dict que tienen un valor de 0.
outliers_dict = {key:val for key, val in outliers_dict.items() if val != 0}

# Graficando el diccionario outliers_dict.
fig, ax = plt.subplots(figsize=(5, 5))
sns.barplot(x=list(outliers_dict.keys()), y=list(outliers_dict.values()))
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
plt.title("Proteinas ribosomales que estan relacionadas a cambios en el ratio 30S/50S")
plt.xlabel("Proteinas ribosomales")
plt.ylabel("Número de condiciones")
plt.savefig(cwd+"/imagenes"+"/outliers_mori.png", dpi=300, bbox_inches='tight')

# Heatmap del dataframe ribo_mori_set_norm, con solo las condiciones que estan en la lista list_pos_outliers.
# Y con solo las proteinas que estan en el diccionario outliers_dict.
# El heatmap muestra los valores de las condiciones que estan por encima o por debajo de la media más o menos 1 desviación estándar.
# El heatmap muestra los valores de las proteinas que estan por encima o por debajo de la media más o menos 2 desviaciones estándar.
fig, ax = plt.subplots(figsize=(5, 5))
sns.heatmap(ribo_mori_set_norm.loc[:,list_pos_outliers].loc[list(outliers_dict.keys())], cmap="coolwarm", ax=ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
plt.title("Proteinas ribosomales que estan relacionadas a cambios en el ratio 30S/50S")
plt.xlabel("Condiciones de crecimiento")
plt.ylabel("Proteinas ribosomales")
plt.savefig(cwd+"/imagenes"+"/heatmap_rp_variaciones.png", dpi=300, bbox_inches='tight')

## Imprimiendo un mensaje de despedida.
print("El programa ha finalizado con éxito.")
print("Las imagenes generadas se encuentran en la carpeta imagenes.")