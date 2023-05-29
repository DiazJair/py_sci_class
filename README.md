# Manual de usuario.

Este manual de usuario contiene las instrucciones basicas para usar las funciones para la obtención de los datos proteomicos obtenidos por Mori et al., 2021.
El artículo se puede consultar en el siguiente link: https://www.embopress.org/doi/full/10.15252/msb.20209536

# Dependencias requeridas:
En la carpeta se encuentra un archivo titulado "dependencias.txt", este archivo contiene las dependencias requeridas y puede usarlo para instalarlas en un ambiente de conda ya existente, o para verificar la versión con un ambiente previo.

Adicionalmente se incluye el archivo env_dup.yml, el cual permite crear un ambiente nuevo con el nombre: py_sci que instalará todas las dependencias requeridas. Para ello solo debe usar el siguiente comando: 
conda env create -f env_dup.yml

# Utilización:

Las funciones estan almacenadas en un script titulado prot_func.py, Las funciones estan debidamente documentas y las puede consultar abriendo el archivo con cualquier lector de texto como nano, vim, kate, etc..

Se incluye un ejemplo practico de utilización en el archivo jair_ernesto_diaz_ramirez.py. Para ejecutarlo solo se tiene que ejecutar el siguiente comando:
python3 jair_ernesto_diaz_ramirez.py
El ejemplo consiste en el analisís de la estequiometria de las proteinas ribosomales. Las imagenes se almacenan en la carpeta imagenes que es creada automaticamente al ejecutar el código. 

# Cualquier duda o comentrario al correo: diaz.jair.e@gmail.com
