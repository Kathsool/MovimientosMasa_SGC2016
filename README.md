# MovimientosMasa_SGC2016
Este código realizado para ArcGIS obtiene la amenaza por movimientos en masa a escala detallada (1:2000/1:5000), según la Guía metodológica para estudios de amenaza, vulnerabilidad y riesgo por movimientos en masa (Avila et al., 2016)


Se tienen 3 archivos:
- Campos.png
- MM_Urbano.tbx
- MM_Urbano_TB.py

El .tbx es la toolbox para abrir en ArcGIS, a la cual se debe asociar el archivo .py (Al abrir la Toolbox le dan en propiedades a MMUrbano, luego en Source y seleccionan o pegan la ruta del .py

El código necesita de un shape de UGIs (Unidades Geologicas para Ingeniería) y el DEM o DTM. Además, solicita otros datos que se encuentran descritos en la Guía del Servicio, como el Coeficiente de aceleración horizontal, y los límites del factor de seguridad. 
Adicionalmente, tiene como opciones para elegir (i) si no se la amenaza en pendientes menores a 5 grados (En caso de que se chulee, las zonas con pendientes menores quedan en amenaza baja); (ii) Si se suaviza con el proceso de ArcGIS de Majority Filter; (iii) y si se tienen shapes de eventos por movimientos en masa, los cuales en caso de que existan se clasifican como amenaza alta.

La imagen Campos.png, tiene los nombres que deben tener las columnas del shape UGIs para que funcione el código.

Suerte! 
Cualquier cosa me escriben nksolera@unal.edu.co
