# Generated by Django 3.1.3 on 2020-11-30 07:51

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('variants', '0003_auto_20201130_0615'),
    ]

    operations = [
        migrations.RenameField(
            model_name='target',
            old_name='basenji_info',
            new_name='info',
        ),
        migrations.RemoveField(
            model_name='target',
            name='cell_type',
        ),
        migrations.RemoveField(
            model_name='target',
            name='epi',
        ),
        migrations.AddField(
            model_name='target',
            name='num',
            field=models.PositiveIntegerField(default=1),
        ),
    ]
