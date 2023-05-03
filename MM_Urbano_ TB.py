import arcpy
from arcpy import env
from arcpy.sa import *
import numpy as np
import math


arcpy.CheckOutExtension("Spatial")
#################################
# # Check out Spacial Analyst
# try:
#     if arcpy.CheckExtension("Spatial") == "Available":
#         arcpy.CheckOutExtension("Spatial")
#         print ("Checked out \"Spatial\" Extension")
#     else:
#         raise LicenseError
# except LicenseError:
#     print ("Spatial Analyst license is unavailable")
# except:
#     print (arcpy.GetMessages(2))

## Código desarrollado por Natalia Katherine Soler Aragon
## 17 Marzo 2023
## Correo: nksolera@unal.edu.co
##############################################

#Carpeta de inicio
origen = arcpy.GetParameterAsText(0)
#Crea carpeta de los intermedios
arcpy.management.CreateFolder(origen, r'01_Apoyo')
pathbas = origen + r'\01_Apoyo'
arcpy.env.workspace = origen

# variables con Paths
ugis = arcpy.GetParameterAsText(1)
mdt = arcpy.GetParameterAsText(2)

## Datos para el Factor Sismico solo NRS-10
Aa = float(arcpy.GetParameterAsText(3))

## Datos para el factor de seguridad
fs1a = arcpy.GetParameterAsText(4)
fs2a = arcpy.GetParameterAsText(5)
fs1 = float(fs1a)
fs2 = float(fs2a)
print(arcpy.AddMessage('Factor de seguridad 1: '+ str(fs1)))
print(arcpy.AddMessage('Factor de seguridad 2: '+ str(fs2)))

bo_p5 = arcpy.GetParameterAsText(6)
bo_sm = arcpy.GetParameterAsText(7)
bo_mm = arcpy.GetParameterAsText(8)
mm = arcpy.GetParameterAsText(9)
bo_sud = arcpy.GetParameterAsText(10)

# Nos deja reescribir archivos
arcpy.env.overwriteOutput = True

### GDBs y paths
arcpy.CreateFileGDB_management(origen, "Insumos.gdb")
arcpy.CreateFileGDB_management(origen, "Intermedios.gdb")
arcpy.CreateFileGDB_management(origen, "Productos.gdb")

insgdb = origen + '\Insumos.gdb'
intgdb = origen + '\Intermedios.gdb'
progdb = origen + '\Productos.gdb'

print(arcpy.AddMessage('Se crearon las GDBs'))

##INTERMEDIOS
fr = intgdb + '\Friccion'
de = intgdb + '\Densidad'
co = intgdb + '\Cohesion'
ta = intgdb + '\TablaAgua'
es = intgdb + '\Espesor'
reag = intgdb + '\RelacionAgua'
coeamp = intgdb + '\CoeficienteAmplificacion'
coeimp = intgdb + '\CoeficienteImportancia'
pe = intgdb + '\Pendiente_grados'
pera = intgdb+ '\Pendiente_radianes'
fsis = intgdb + '\FactorSismico'
pen = intgdb + '\Pendiente_IGAC'
Pen_me_5 = intgdb + '\Pendiente_Menor5'

##BAS
#Coeficiente de ampliacion sismica COAS
coas = pathbas + '\St_COAs.tif'
coas_shd = pathbas + '\St_CoAsShD.shp'
stf = pathbas + '\STFinal.tif'
cosc = pathbas + '\CosCos.tif'
sicos = pathbas + '\SinCos.tif'
fse_rec = pathbas + '\FSe_reclass.tif'
fse_rec_sin = pathbas + '\FSe_reclass_sinSmooth.tif'
fse_shp = pathbas + '\FSe_reclass.shp'
pe_rec = pathbas + '\Pend5_reclass.tif'
pe5r_sh = pathbas + '\Pend5_reclass.shp'
pe5_sh = pathbas + '\Pendiente_Menor5.shp'
fse_pe5_er = pathbas + '\FSeg_EraseMenor5.shp'
fse_mmer = pathbas + '\FSeg_Erase_MM.shp'
cr_fs_sus_r = pathbas + '\FSe_Sus_reclass.tif'
cr_fs_sis_r = pathbas + '\FSe_Sis_reclass.tif'
cr_fs_ta_r = pathbas + '\FSe_Tab_reclass.tif'
ugisclip = pathbas + '\Ugis_buffer.shp'
# Detonantes y Suscep
cr_fs_sus_rs= pathbas + '\FSe_sus_reclass.shp'
cr_fs_sis_rs= pathbas + '\FSe_sis_reclass.shp'
cr_fs_ta_rs= pathbas + '\FSe_ta_reclass.shp'


####### PRODUCTOS
facsegu = progdb + '\FactorSeguridad'
fse_rdi = progdb + '\Factor_Seguridad'
fse_pe5sc = progdb + '\Factor_Seguridad_P5_sinclip'
fse_pe5 = progdb + '\Factor_Seguridad_P5'
fse_mm = progdb + '\Factor_Seguridad_MM'
fse_pe5_mm = progdb + '\Factor_Seguridad_P5_MM'
# susceptibilidad y detonantes
cr_fs_sus = progdb + '\FactorSeguridad_Suscep'
cr_fs_sis = progdb + '\FactorSeguridad_Sismico'
cr_fs_ta = progdb + '\FactorSeguridad_Tabla'
cr_fs_sus_di = progdb + '\Factor_Seguridad_Susceptib'
cr_fs_sis_di = progdb + '\Factor_Seguridad_Sismo'
cr_fs_ta_di = progdb + '\Factor_Seguridad_TablaAgua'


##INSUMOS
ugisgdb = insgdb + '\UgisGDB'
mdtcl = insgdb + '\Mdt_clip'
arcpy.management.CopyFeatures(ugis, ugisgdb)
arcpy.analysis.Buffer(ugisgdb, ugisclip, "15 Meters")

maskmdt = ExtractByMask(mdt, ugisclip)
maskmdt.save(mdtcl)

arcpy.conversion.RasterToGeodatabase([Raster(mdt)], insgdb)



###  Relacion tabla de agua
arcpy.management.AddField(ugisgdb, 'RelacionAgua', 'DOUBLE')
# Llenar tabla de atributos
with arcpy.da.UpdateCursor(ugisgdb, ['Espesor','TablaAgua', 'RelacionAgua']) as tablecursor:
    for row in tablecursor:
        row[2] = row[1]/row[0]
        tablecursor.updateRow(row)

print(arcpy.AddMessage('Se calculo la relacion de la tabla de agua en el shape.'))

# Tamano celda
cellsizex = arcpy.GetRasterProperties_management(mdtcl, 'CELLSIZEX')
cellsizey = arcpy.GetRasterProperties_management(mdtcl, 'CELLSIZEY')

# Raster to polygon
arcpy.conversion.PolygonToRaster(ugisgdb, 'Friccion', fr, 'CELL_CENTER', 'NONE', cellsizex)
arcpy.conversion.PolygonToRaster(ugisgdb, 'Densidad', de, 'CELL_CENTER', 'NONE', cellsizex)
arcpy.conversion.PolygonToRaster(ugisgdb, 'Cohesion', co, 'CELL_CENTER', 'NONE', cellsizex)
arcpy.conversion.PolygonToRaster(ugisgdb, 'TablaAgua', ta, 'CELL_CENTER', 'NONE', cellsizex)
arcpy.conversion.PolygonToRaster(ugisgdb, 'Espesor', es, 'CELL_CENTER', 'NONE', cellsizex)
arcpy.conversion.PolygonToRaster(ugisgdb, 'RelacionAgua', reag, 'CELL_CENTER', 'NONE', cellsizex)
arcpy.conversion.PolygonToRaster(ugisgdb, 'CAmplif', coeamp, 'CELL_CENTER', 'NONE', cellsizex)
arcpy.conversion.PolygonToRaster(ugisgdb, 'CImport', coeimp, 'CELL_CENTER', 'NONE', cellsizex)
####################################################

### Calculo de la pendiente según IGAC
# Pendiente en grados
outSlope = Slope(mdtcl, "DEGREE")
#outSlope.save(pera)
max_pe = arcpy.GetRasterProperties_management(outSlope, 'MAXIMUM')
min_pe = arcpy.GetRasterProperties_management(outSlope, 'MINIMUM')
arcpy.AddMessage("El valor maximo de pendiente es: " + str(max_pe))
arcpy.AddMessage("El valor minimo de pendiente es: " + str(min_pe))

pe_recl = Reclassify(outSlope, 'Value',
                       RemapRange([[0,2,1], [2, 4, 2],[4, 7, 3], [7, 14, 4],
                                   [14, 27, 5], [27, 37, 6],[37, 100, 7],
                                   ]))
#pe_recl.save(pe_re)
arcpy.conversion.RasterToPolygon(pe_recl, pen, 'SIMPLIFY', 'Value', 'MULTIPLE_OUTER_PART')
# En la tabla de Atributos
arcpy.management.AddField(pen, 'Pendiente', 'TEXT')
arcpy.management.AddField(pen, 'Rango', 'TEXT')

maxi_pe = round(float(str(max_pe)),2)
mini_pe = round(float(str(min_pe)),2)
pen_values = [str(mini_pe)+' a 2', '2 a 4', '4 a 7', '7 a 14', '14 a 27', '27 a 37','37 a '+str(maxi_pe) ]
pen_rango = ['Plano', 'Ligeramente inclinado', 'Moderadamente inclinado', 
             'Fuertemente inclinado', 'Ligeramente escarpado', 'Moderadamente escarpado',
               'Fuertemente escarpado']
yy = 0
# Llenar tabla de atributos
with arcpy.da.UpdateCursor(pen, ['Pendiente', 'Rango']) as cursor:
    for row in cursor:
        row[0] = pen_values[yy]
        row[1] = pen_rango[yy]
        yy = yy +1
        cursor.updateRow(row)


####################################################
# Pendiente en grados
outSlope = Slope(mdtcl, "DEGREE")
outSlope.save(pe)
# Pendiente a radianes
radslope2 = outSlope * math.pi/180
radslope2.save(pera)


## ST (Factor sismico > Coeficiente de ampliacion sismica) - coas
max_slope = arcpy.GetRasterProperties_management(outSlope, 'MAXIMUM')
#print(arcpy.AddMessage('La pendiente maxima es: '+str(max_slope)))
coas_calculo = Reclassify(outSlope, 'Value', RemapRange([[0,15,1],[15,30,2],[30,100,3]]))
coas_calculo.save(coas)
arcpy.conversion.RasterToPolygon(coas_calculo, coas_shd, 'SIMPLIFY', 'Value','MULTIPLE_OUTER_PART')

# En la tabla de Atributos
arcpy.management.AddField(coas_shd, 'ST', 'DOUBLE')
st_values = [1.0, 1.2, 1.4]
y = 0
# Llenar tabla de atributos
with arcpy.da.UpdateCursor(coas_shd, 'ST') as cursor:
    for row in cursor:
        row[0] = st_values[y]
        y = y +1
        cursor.updateRow(row)
# ST Final a shape
arcpy.conversion.PolygonToRaster(coas_shd, 'ST', stf, 'CELL_CENTER', 'NONE', cellsizex)
print(arcpy.AddMessage('Se guardo stf (Coeficiente de ampliacion sismica)'))


# Calculo  FACTOR SISMICO
cal_sis = Raster(stf) * Aa * Raster(coeamp) * Raster(coeimp)
cal_sis.save(fsis)

print(arcpy.AddMessage('Se calculo el Factor Sismico'))

## Necesarios para el factor de seguridad
coscos = Cos(pera)*Cos(pera)
coscos.save(cosc)
sincos = Sin(pera)*Cos(pera)
sincos.save(sicos)
tang = Tan(fr)

## a raster
rco = Raster(co)
rfr = Raster(fr)
rde = Raster(de)
rfs = Raster(fsis)
esp = Raster(es)
rag = Raster(reag)

facseg = (rco +((rde * esp * coscos)-(rfs * esp * rde * sincos)-(rag * esp * coscos))*tang)/((rde*esp*sincos)+(rfs*rde*esp*coscos))
facseg.save(facsegu)

print(arcpy.AddMessage('Se calculo el Factor de Seguridad'))

#Limites del Factor sismico
max_facseg = arcpy.GetRasterProperties_management(facseg, 'MAXIMUM')
print(arcpy.AddMessage('El maximo del Factor de Seguridad es: '+ str(max_facseg)))
min_facseg = arcpy.GetRasterProperties_management(facseg, 'MINIMUM')
print(arcpy.AddMessage('El minimo del Factor de Seguridad es: '+ str(min_facseg)))
facseg_rec = Reclassify(facseg, 'Value', RemapRange([[0, fs1a, 1], [fs1a, fs2a, 2], [fs2a, 800000, 3]]))
f_max = max_facseg.getOutput(0)

if float(f_max) >= 800000:
    print(arcpy.AddMessage('ERROR: El factor de seguridad es mayor al limite maximo establecido, aunque salgan resultados son erroneos ya que tendran huecos!'))
else: pass

#####################################################################################
######### CON Smooth
if bo_sm == 'true':
    print(arcpy.AddMessage('Se eligio hacer Majority Filter, se esta procesando'))
    facseg_rec.save(fse_rec_sin)
    majo = MajorityFilter(fse_rec_sin, 'EIGHT', 'MAJORITY')
    majo.save(fse_rec)
elif bo_sm == 'false':
    # Sin smooth
    print(arcpy.AddMessage('No se eligio hacer Smooth'))
    facseg_rec.save(fse_rec)
############################################################################
#####################
## Raster a shape FACTOR SEGURIDAD
arcpy.conversion.RasterToPolygon(fse_rec, fse_shp, 'SIMPLIFY', 'Value')
arcpy.management.Dissolve(fse_shp, fse_rdi, 'gridcode')
arcpy.management.AddField(fse_rdi, 'FacSegur', 'TEXT')
arcpy.management.AddField(fse_rdi, 'Amenaza', 'TEXT')
maxi = round(float(str(max_facseg)),2)
mini = round(float(str(min_facseg)),2)
facse_values = ['0 a '+ str(fs1),str(fs1)+' a '+str(fs2),'Mayor a '+ str(fs2)]
ame_values = ['Alta', 'Media', 'Baja']

# Llenar tabla de atributos
with arcpy.da.UpdateCursor(fse_rdi, ['FacSegur', 'Amenaza','gridcode']) as cursor:
    for row in cursor:
        existing_value = row[2]

        # Amenaza Alta
        if existing_value == long(1):
            row[0] = facse_values[0]
            row[1] = ame_values[0]
            cursor.updateRow(row)
        # Amenaza Media
        if existing_value == long(2):
            row[0] = facse_values[1]
            row[1] = ame_values[1]
            cursor.updateRow(row)
        # Amenaza Baja
        if existing_value == long(3):
            row[0] = facse_values[2]
            row[1] = ame_values[2]
            cursor.updateRow(row)


###################################################################
###################################################################
###################################################################
# Susceptibilidad y detonantes:
if bo_sud == 'true':
    print(arcpy.AddMessage('Se eligio calcular la Susceptibilidad y los detonantes'))
    rfs1 = 0
    rag1 = 0
    # Susceptibilidad
    c_fs_sus = (rco +((rde * esp * coscos)-(rfs1 * esp * rde * sincos)-(rag1 * esp * coscos))*tang)/((rde*esp*sincos)+(rfs1*rde*esp*coscos))
    c_fs_sus.save(cr_fs_sus)
    # sismico
    c_fs_sis = (rco +((rde * esp * coscos)-(rfs * esp * rde * sincos)-(rag1 * esp * coscos))*tang)/((rde*esp*sincos)+(rfs*rde*esp*coscos))
    c_fs_sis.save(cr_fs_sis)
    # tabla de agua
    c_fs_ta = (rco +((rde * esp * coscos)-(rfs1 * esp * rde * sincos)-(rag * esp * coscos))*tang)/((rde*esp*sincos)+(rfs1*rde*esp*coscos))
    c_fs_ta.save(cr_fs_ta)
else: pass

# Susceptibilidad y detonantes:
sustext = ['Susceptibilidad', 'Sismico', 'TablaAgua']
susdet = [cr_fs_sus, cr_fs_sis, cr_fs_ta]
sd_fs_rec = [cr_fs_sus_r, cr_fs_sis_r, cr_fs_ta_r]
sd_fs_shp = [cr_fs_sus_rs, cr_fs_sis_rs, cr_fs_ta_rs]
sd_fs_dis = [cr_fs_sus_di, cr_fs_sis_di, cr_fs_ta_di]
if bo_sud == 'true':
    for x_sd in range(len(susdet)):
        #Limites del Factor de Seguridad
        max_fs_sd = arcpy.GetRasterProperties_management(susdet[x_sd], 'MAXIMUM')
        min_fs_sd = arcpy.GetRasterProperties_management(susdet[x_sd], 'MINIMUM')
        f_maxf = max_fs_sd.getOutput(0)
        print(arcpy.AddMessage('Procesando: {}'.format(sustext[x_sd])))
        if float(f_maxf) >= float(800000):
            print(arcpy.AddMessage('ERROR: El factor de seguridad es mayor al limite maximo establecido, aunque salgan resultados son erroneos ya que tendran huecos!'))
        else: pass

        # Reclassify del Factor de Seguridad
        sd_fs_re = Reclassify(susdet[x_sd], 'Value', RemapRange([[0, fs1a, 1], [fs1a, fs2a, 2], [fs2a, 800000, 3]]))
        sd_fs_re.save(sd_fs_rec[x_sd])
        
        ## Raster a shape FACTOR SEGURIDAD
        arcpy.conversion.RasterToPolygon(sd_fs_rec[x_sd], sd_fs_shp[x_sd], 'SIMPLIFY', 'Value')
        arcpy.management.Dissolve(sd_fs_shp[x_sd], sd_fs_dis[x_sd], 'gridcode')
        arcpy.management.AddField(sd_fs_dis[x_sd], 'FacSegur', 'TEXT')
        arcpy.management.AddField(sd_fs_dis[x_sd], 'Susceptibilidad', 'TEXT')
        maxi = round(float(str(max_fs_sd)),2)
        mini = round(float(str(min_fs_sd)),2)
        facse_values = ['0 a '+ str(fs1),str(fs1)+' a '+str(fs2), 'Mayor a '+ str(fs2)]
        ame_values = ['Alta', 'Media', 'Baja']
        

        # Llenar tabla de atributos
        with arcpy.da.UpdateCursor(sd_fs_dis[x_sd], ['FacSegur', 'Susceptibilidad','gridcode']) as cursor:
            for row in cursor:
                existing_value = row[2]

                # Amenaza Alta
                if existing_value == long(1):
                    row[0] = facse_values[0]
                    row[1] = ame_values[0]
                    cursor.updateRow(row)
                # Amenaza Media
                if existing_value == long(2):
                    row[0] = facse_values[1]
                    row[1] = ame_values[1]
                    cursor.updateRow(row)
                # Amenaza Baja
                if existing_value == long(3):
                    row[0] = facse_values[2]
                    row[1] = ame_values[2]
                    cursor.updateRow(row)
elif bo_sud == 'false':
    print(arcpy.AddMessage('No se eligio calcular Susceptibilidad y detonantes'))
###################################################################
###################################################################
###################################################################



######### Pendientes menores a  5 grados
if bo_p5 == 'true':
    print(arcpy.AddMessage('Se eligio quitar las pendientes menores a 5 grados'))
    pe_re = Reclassify(pe, 'Value', RemapRange([[0,5,1],[5,100,2]]))
    pe_re.save(pe_rec)
    arcpy.conversion.RasterToPolygon(pe_rec, pe5r_sh, 'SIMPLIFY', 'Value')
    arcpy.management.Dissolve(pe5r_sh, pe5_sh, 'gridcode')
    arcpy.management.AddField(pe5_sh, 'Pendiente', 'TEXT')
    arcpy.management.AddField(pe5_sh, 'Amenaza', 'TEXT')
    arcpy.management.AddField(pe5_sh, 'FacSegur', 'TEXT')
    pend5_values = ['Pendiente menor a 5', 'Mayor a 5']
    ame_values = ['Baja', 'NN']
    y = 0
    # Llenar tabla de atributos
    with arcpy.da.UpdateCursor(pe5_sh, ['Pendiente', 'Amenaza', 'FacSegur']) as cursor:
        for row in cursor:
            row[0] = pend5_values[y]
            row[1] = ame_values[y]
            row[2] = pend5_values[y]
            y = y +1
            cursor.updateRow(row)
    #Shape solo menores a 5
    arcpy.MakeFeatureLayer_management(pe5_sh,'pe5_layer')
    arcpy.management.SelectLayerByAttribute('pe5_layer', 'NEW_SELECTION',""" "Pendiente" = 'Pendiente menor a 5' """)
    arcpy.conversion.FeatureClassToFeatureClass('pe5_layer', intgdb, 'Pendiente_Menor5')
    # Erase and merge
    arcpy.analysis.Erase(fse_rdi, Pen_me_5, fse_pe5_er)
    arcpy.management.Merge([fse_pe5_er,Pen_me_5], fse_pe5sc)
    arcpy.analysis.Clip(fse_pe5sc, ugisgdb, fse_pe5)
    arcpy.management.Delete(fse_pe5sc)


elif bo_p5 == 'false':
    print(arcpy.AddMessage('No se eligio quitar las pendientes menores a 5 grados'))

###############################################
## Movimientos en masa 
if bo_mm == 'true':
    print(arcpy.AddMessage('Si hay movimientos en masa'))
    mmgdb = insgdb + '\MMGDB'
    mmgdbd = pathbas + '\MM_diss.shp'
    arcpy.management.CopyFeatures(mm, mmgdb)
    
    # MM
    arcpy.management.AddField(mmgdb, 'FacSegur', 'TEXT')
    arcpy.management.AddField(mmgdb, 'Amenaza', 'TEXT')
    # Llenar tabla de atributos
    with arcpy.da.UpdateCursor(mmgdb, ['FacSegur', 'Amenaza']) as cursor:
        for row in cursor:
            row[0] = 'MM reportados'
            row[1] = 'Alta'
            cursor.updateRow(row)
    arcpy.management.Dissolve(mmgdb, mmgdbd, ['FacSegur', 'Amenaza'])

    ### ERASE
    # Sin pendientes 
    if bo_p5 == 'false':
        arcpy.analysis.Erase(fse_rdi, mmgdbd, fse_mmer)                     
        arcpy.management.Merge([mmgdbd, fse_mmer], fse_mm)
        
    # Con pendientes
    elif bo_p5 == 'true':
        arcpy.analysis.Erase(fse_pe5, mmgdbd, fse_mmer)                    
        arcpy.management.Merge([mmgdbd, fse_mmer], fse_pe5_mm)

elif bo_mm == 'false':
    print(arcpy.AddMessage('No hay movimientos en masa :)'))

print(arcpy.AddMessage('Acabamos :)'))

## Código desarrollado por Natalia Katherine Soler Aragon
## 17 Marzo 2023
## Correo: nksolera@unal.edu.co