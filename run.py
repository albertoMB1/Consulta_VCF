# LUIS ALBERTO MARTÍNEZ BRAVO
import os
import pandas as pd
import streamlit as st

import pandas as pd

from io import StringIO
import json
from LIBRERIA import vcf41 as vc



# Configuración página inicial
st.set_page_config(
                    page_title="TFG UPV",
                    initial_sidebar_state="expanded",
                    layout='wide'
                    )
# Paginas de la app
PAGES = [
    'ANALISIS DATOS',
    'BUSQUEDA',
    'ARCHIVO'
]

# Diccionario radio button
dic ={'1': 'ANALISIS DATOS',
          '2': 'BUSQUEDA',
          '3': 'ARCHIVO' }


RUTA = 'VCF/GRCh38_latest_clinvar.vcf'


# Carpeta Base
ruta_base_1 = 'Datos'
# Lista de subcarpetas para crear dentro de la carpeta 'Datos'
subcarpetas = ['Cromosoma', 'CLNDNINCL', 'CLNDN','SUBCAM']

# Crear la carpeta 'Datos' si no existe
os.makedirs(ruta_base_1, exist_ok=True)

# Crear las subcarpetas dentro de 'Datos' si no existen
for carpeta in subcarpetas:
    os.makedirs(os.path.join(ruta_base_1, carpeta), exist_ok=True)





# Carga el diccionario de clndn para encontrar más rápido las enfermedades
# Es un diccionario cuya clave son listas que contiene las líneas del archivo VCF
# Donde se encuentra

archivo_indices_cldnd = 'DATOS/CLNDN/CLNDN.JSON'

if os.path.exists(archivo_indices_cldnd):
    with open(archivo_indices_cldnd, 'r') as fh:
        indices_enfermedades_clndn = json.load(fh)
else:
    st.error("ATENCIÓN: Ejecuta PRIMERO análisis -> CLNDN ")

archivo_indices_clndnincl = 'DATOS/CLNDNINCL/CLNDNINCL.JSON'

if os.path.exists(archivo_indices_clndnincl):
    with open(archivo_indices_clndnincl, 'r') as fh:
        indices_enfermedades_clndnincl = json.load(fh)
else:
    st.error("ATENCIÓN: Ejecuta PRIMERO análisis -> CLNDNINCL ")

if not os.path.exists('VCF/GRCh38_latest_clinvar.vcf'):
    st.error(("ATENCIÓN: No está el archivo VCF con el nombre adecuado"))



# Define las columnas específicas
columns = ["#CHROM", "#POS", "#ID", "#REF", "#ALT", "#QUAL", "#FILTER", "#INFO"]



def run():

    st.sidebar.title('Luis Alberto MB')
    if st.session_state.page:
        page=st.sidebar.radio('## Navegación', PAGES, index=st.session_state.page)
    else:
        page=st.sidebar.radio('## Navegación', PAGES, index=1,)
    

    if page == 'ANALISIS DATOS':
        st.title("ANÁLISIS Y CREACION DE DICCIONARIOS PARA LOS OTROS APARTADOS")
        st.sidebar.write("""
            ## Sobre:
            Acciones que se ejecutan:\n
                    """)
        st.sidebar.write("""- Análisis de ID Único en VCF
                         \n- Extracción valores de enfermedades CLNDN, CLNVI
                        \n- Valores de CLNDN
                         \n- Número de Cromosomas  
                         """)
        subcolumna1, subcolumna2 = st.columns(2)
        with subcolumna1:
            st.session_state.data_type = st.radio("Elige una opción", ('ID', 'Cromosoma', 'CLNDN','CLNDNINCL', 'Subcampo_de_Info'), index = 0)
        
        if st.session_state.data_type =='ID':

            if st. button("Verificar ID Distintos"):
                numero1 = vc.numero_datos_vcf(RUTA)
                numero2 = len(vc.obtain_all_ID_inVCF(RUTA))

                if numero1 == numero2:
                    st.success(f"Enhorabuena. Los ID son únicos\nNumero datos = {numero1}\nNúmero registros distintos = {numero2}")
                else:
                    st.error("Los ID no son únicos")

        
    
        if st.session_state.data_type =='Cromosoma':
                
            if st.button("Chromosomas"):
                
                data = sorted(list(vc.obtain_all_names_chrom_inVCF(RUTA)))

                df = pd.DataFrame([data], columns=None)
                st.dataframe(df)

                cont = vc.numero_datos_vcf(RUTA)
                st.success(f"Número de datos = {cont} .")

                dic = {clave: 0 for clave in data}
                ruta_dic_c ='DATOS/Cromosoma/CROMOSOMAS.JSON'
                
                crea_directori_si_no_existe(os.path.dirname(ruta_dic_c))

                vc.guardar_dicc_json(ruta_dic_c, dic)


        if st.session_state.data_type == 'CLNDN':
            if st.button('Busca'):

                #data = sorted(list(vc.busca_todos_valores_en_un_subcampo_de_info(ruta,'CLNDN')))
                dicc = {}
                with open(RUTA, 'r') as arch:
                    num = 0
                    for linea in arch:
                        num +=1
                        if not linea.startswith('#'):
                            columnas = linea.split('\t')
                            info_dict = parse_info_column(columnas[7])
                            
                            clndn = info_dict.get('CLNDN')
                            if clndn:
                                update_dictionary(dicc, clndn, num)

                data = sorted(list(dicc.keys()))
               


                ruta_dic ='DATOS/CLNDN/CLNDN.JSON'
                crea_directori_si_no_existe(os.path.dirname(ruta_dic))

                vc.guardar_dicc_json(ruta_dic, dicc)
                st.success(f'Guardado en disco en {ruta_dic}\nNúmero de enfermedades: {len(dicc)}',icon='✅')
                df = pd.DataFrame(data, columns=["Enfermedad"])
                st.dataframe(df)

        if st.session_state.data_type == 'Subcampo_de_Info':
            if st.button('Busca'):
                data = sorted(list(vc.encuentra_param_vcf4(RUTA)))
                dic = {clave: 0 for clave in data}

                ruta_dic ='DATOS/SUBCAM/SUBCAM.JSON'
                crea_directori_si_no_existe(os.path.dirname(ruta_dic))
 
                vc.guardar_dicc_json(ruta_dic, dic)
                st.success(f'Guardado en disco en {ruta_dic}\nNúmero de subcampos: {len(dic)}')
                df = pd.DataFrame(data, columns=["SUBCAMPOS"])
                st.dataframe(df)


        if st.session_state.data_type == 'CLNDNINCL':
            if st.button('Busca'):

                dicc = {}
                with open(RUTA, 'r') as arch:
                    num = 0
                    for linea in arch:
                        num +=1
                        if not linea.startswith('#'):
                            columnas = linea.split('\t')
                            info_dict = parse_info_column(columnas[7])
                            
                            clndn = info_dict.get('CLNDNINCL')
                            if clndn:
                                update_dictionary(dicc, clndn, num)

                data = sorted(list(dicc.keys()))


                vc.guardar_dicc_json('DATOS/CLNDNINCL/CLNDNINCL.JSON', dicc)
                st.success(f'Guardado en disco en DATOS/CLNDNINCL/CLNDNINCL.JSON\nNúmero de enfermedades: {len(dicc)}')
                df = pd.DataFrame(data, columns=["Enfermedad"])
                st.dataframe(df)


    if page == 'BUSQUEDA':
        st.sidebar.write("""
            ## Sobre:
            Búsqueda de diferentes campos:\n
                    """)
        st.sidebar.write("""
                         \n- Introduciendo ID
                        \n- Por enfermedad en CLNDN
                         \n- Pos
                         """)
        st.title("Búsquedas")

        subcolom1, subcolom2 = st.columns(2)
        with subcolom1:
            st.session_state.data_type = st.radio("Elige una ítem para buscar por:", ('ID',  'POS', 'Enfermedad_CLNDN', 'Enfermedad_CLNDNINCL'), index = 0)

        # Búsqueda por CLNDN
        if st.session_state.data_type =='Enfermedad_CLNDN':
            
            ruta_archivo = "DATOS/CLNDN/CLNDN.json"
            diccionario_leido = vc.leer_dic_json(ruta_archivo)
            opcion_seleccionada = st.selectbox("Selecciona una enfermedad", list(diccionario_leido.keys()))
            if st.button("Busca por CLNDN"):
            # Convertir result.stdout a DataFrame
                _, data, _ = mostrar_dataframe(opcion_seleccionada,'CLNDN')
                st.download_button(label='Descargar',data=data.getvalue().encode(),file_name=f'Resultados_CLNDN_{elimina_barra(opcion_seleccionada)}.txt')

        # Búsqueda por ID
        if st.session_state.data_type =='ID':
            entrada = st.text_input("Escribe número ID")
            if st.button("Buscando por ID"):
                mostrar_dataframe_id(entrada)
        # Búsqueda por POS
        if st.session_state.data_type =='POS':
            entrada = st.text_input("Escribe número POS")
            if st.button("Bucar posición"):
                mostrar_dataframe_pos(entrada)
                    
        # Búsqued por CLNDNINCL
        if st.session_state.data_type =='Enfermedad_CLNDNINCL':
                dic_CLNDNINCL = vc.leer_dic_json('DATOS/CLNDNINCL/CLNDNINCL.json')
                opcion_seleccionada = st.selectbox("Selecciona enfermedad", list(dic_CLNDNINCL))
                if st.button("Buscar"):
                    _ , data, _ = mostrar_dataframe(opcion_seleccionada,'CLNDNINCL')
                    st.download_button(label='Descargar',data=data.getvalue().encode(),file_name=f'Resultados_CLNDNINCL_{elimina_barra(opcion_seleccionada)}.txt')
     
                    
    if page == 'ARCHIVO':
            st.sidebar.write("""
            ## Sobre:
            Creación del archivo VCF por patología\n
                    """)
            #st.title("Create VCF FILE")
            st.subheader("Crea archivo VCF filtrado según CLNDN")
            ruta_archivo = "DATOS/CLNDN/CLNDN.json"
            diccionario_leido = vc.leer_dic_json(ruta_archivo)
            opcion_seleccionada = st.selectbox("Selecciona una enfermedad", list(diccionario_leido.keys()))
            
            if st.button("Buscar"):
            # Convertir result.stdout a DataFrame
            
                _, _, busqueda = mostrar_dataframe(opcion_seleccionada,'CLNDN')
                print(opcion_seleccionada)

                arch_ruta= f'Output/{elimina_barra(opcion_seleccionada)}.vcf'
     
                crear_archivo_vcf(busqueda, arch_ruta)

                st.success(f"Archivo generado correctamente. Nombre = Ouput_{arch_ruta}",icon='✅')
                descargar(arch_ruta)
                #os.remove(arch_ruta)
def parse_info_column(info_col):
    return {key: value for part in info_col.split(';') if '=' in part for key, value in [part.split('=', 1)]}



#@st.cache_data               
def mostrar_dataframe(opcion_seleccionada, param_str:str):

    #aux = vc.buscar_campo_en_vcf_string(ruta, opcion_seleccionada, param_str)

    busqueda=''
    if param_str=='CLNDN':
        busqueda = busca_en_arch(opcion_seleccionada,param_str)

    elif param_str=='CLNDNINCL':
        busqueda = busca_en_arch( opcion_seleccionada, param_str)

    st.toast(f'Búsqueda finalizada', icon='✅')
    st.toast(f'Mostrando Datos',icon='🧑‍💻')
    data = StringIO(busqueda)
    df = pd.read_csv(data, sep="\t", names=columns, header=None, dtype={0: str})  

    st.success(f'Número de resultados = {len(df)}')
    st.dataframe(df)

    return df, data, busqueda

@st.cache_data   
def mostrar_dataframe_pos(opcion_seleccionada:str):
    sol = vc.buscar_pos_en_vcf_str(RUTA, opcion_seleccionada)
    if sol is not None:
        data = StringIO(sol)
        df = pd.read_csv(data, sep="\t", names=columns, header=None, dtype={0: str})
        st.success(f'Número de resultados = {len(df)}')
        st.dataframe(df)
    else:
        st.error(f"No existe la POS = {opcion_seleccionada}")

@st.cache_data  
def mostrar_dataframe_id(opcion_seleccionada:str):
    sol = vc.buscar_id_en_vcf_str(RUTA, opcion_seleccionada)

    if sol is not None:
        data = StringIO(sol)
        df = pd.read_csv(data, sep="\t", names=columns, header=None, dtype={0: str})
        st.success(f'Número de resultados = {len(df)}')
        st.dataframe(df)
    else:
        st.error(f"No existe la ID = {opcion_seleccionada}")


def busca_en_arch(opcion_seleccionada,param_str):
    busqueda=''
    if param_str == 'CLNDN':
        lineas_deseadas = set(indices_enfermedades_clndn[opcion_seleccionada])
    else:
        lineas_deseadas = set(indices_enfermedades_clndnincl[opcion_seleccionada])
    with open(RUTA, 'r') as arch:
        for indice, linea in enumerate(arch,1):
            if indice in lineas_deseadas:
                busqueda += linea
    return busqueda


def crear_archivo_vcf(busqueda, ruta_salida):
    with st.spinner('Guardando archivo... Esto puede tardar unos momentos.'):
        cabecera = vc.extraer_cabeceras_vcf(RUTA)
        cabecera += busqueda
        with open(ruta_salida, 'w') as fd:
            fd.write(cabecera)

    
def descargar(ruta_archivo):
    with open(ruta_archivo, 'r') as f:
        btn = st.download_button(
            label="Descargar archivo VCF",
            data=f,
            file_name=ruta_archivo,
            mime='text/vcf'
        )

def crea_directori_si_no_existe(dir_path):
    """Create directory if it does not exist."""
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def elimina_barra(opcion_seleccionada):
    opcion=''
    for caracter in opcion_seleccionada:
        if caracter in ['/']:
            opcion += '['
        else:
            opcion +=caracter
    return opcion

def parse_info_column(info_col):
    return {key: value for part in info_col.split(';') if '=' in part for key, value in [part.split('=', 1)]}

def update_dictionary(dicc, values, linea_num):
    for value in values.split('|'):
        if value in dicc:
            dicc[value].append(linea_num)
        else:
            dicc[value] = [linea_num]


if __name__ == '__main__':

    if not os.path.exists('Output'):
        os.makedirs('Output')
    if 'loaded' not in st.session_state:
        st.session_state.page = PAGES.index(dic['1'])
        st.session_state['loaded'] = False

    run()