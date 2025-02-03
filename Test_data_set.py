import pandas as pd  
data = pd.read_csv("sird_dataset.csv")
print("Taille des données du fichier CSV :", data.shape)
print("Noms des colonnes :", data.columns)

