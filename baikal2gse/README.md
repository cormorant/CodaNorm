baikal2gse
===

**Программа для конвертации данных из формата Байкал-5 (ХХ) в GSE2**

### Abstract

Waveforms of earthquakes recorded are stored in a regional format Baykal-5. This script converts data from Baykal-5 format to GSE2 format.



### Использование программы

Просмотр всех опций доступен с помощью команды

**baikal2gse.exe -h**



baikal2gse.py [-h] [-V] [-o OUTDIR] dirs [dirs ...]

обязательные аргументы:
  dirs                  путь к данным в формате Байкал

опциональные аргументы:
  -h, --help            справка
  -V, --version         версия программы
  -o OUTDIR, --outdir OUTDIR
                        путь для сохранения выходных данных (по умолчанию "gse2")
