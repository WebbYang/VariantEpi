# VariantEpi

Deploy site: https://variantepi.c4lab.tw


## Setup

conda:
``` bash
conda create --name tf1 python=3.7.4
conda activate tf1
pip install -r requirements.txt
```

docker:
``` bash
docker build . -t c4lab/variantepi
```


## Initialize

``` bash
# create db.sqlite3
docker run -it --rm -v $PWD:/app c4lab/variantepi python manage.py migrate

# download hg19
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
docker run -it --rm -v $PWD:/app -w /app quay.io/biocontainers/samtools:1.15.1--h1170115_0 samtools faidx hg19.fa

# create data folder
touch subitem_itvl.bed
mkdir fig_pred_out model_data
```

If you want to host on your domain, change the hostname in the setting `VariantEpi/settings.py`.
```
ALLOWED_HOSTS = [
    'variantepi.c4lab.tw'
]
```


## Run

``` bash
docker run -it --rm \
    -p 8000:8000 \
    -e DEBUG=false \
    -e SECRET_KEY=a_random_and_secret_string_you_need_to_set \
    -v $PWD/hg19.fa:/app/data/hg19.fa:ro \
    -v $PWD/hg19.fa.fai:/app/data/hg19.fa.fai:ro \
    -v $PWD/db.sqlite3:/app/db.sqlite3 \
    -v $PWD/subitem_itvl.bed:/app/data/subitem_itvl.bed \
    -v $PWD/fig_pred_out:/app/variants/static/pred_out \
    -v $PWD/model_data:/root/.kipoi \
    c4lab/variantepi
```
