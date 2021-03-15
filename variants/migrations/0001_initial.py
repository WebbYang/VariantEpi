# Generated by Django 3.1.3 on 2020-11-30 05:52

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Snp',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('rsid', models.CharField(max_length=32)),
                ('chrm', models.CharField(max_length=32)),
                ('pos', models.PositiveIntegerField()),
                ('ref', models.CharField(max_length=1)),
                ('alt', models.CharField(max_length=1)),
            ],
        ),
        migrations.CreateModel(
            name='target',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('basenji_info', models.CharField(max_length=64)),
                ('cell_type', models.CharField(max_length=32)),
                ('epi', models.CharField(max_length=32)),
            ],
        ),
    ]
