from django.db import models

# Create your models here.
class Snp(models.Model):
    id = models.AutoField(primary_key=True)
    rsid = models.CharField(max_length=32)
    '''chrm = models.CharField(max_length=32,blank=True)
    pos = models.PositiveIntegerField(blank=True)
    ref = models.CharField(max_length=1,blank=True)
    alt = models.CharField(max_length=1,blank=True)'''

    def __str__(self):
        return f"{self.rsid}" # {self.chrm}:{self.pos}

class target(models.Model):
    id = models.AutoField(primary_key=True)
    info = models.CharField(max_length=64)
    num = models.PositiveIntegerField(default=1)
    cell = models.CharField(max_length=64)
    #epi = models.CharField(max_length=32)

    def __str__(self):
        return f"{self.info}: {self.num} in {self.cell}"

#class Score(models.Model):
    #id = models.ForeignKey(Snp)
    #rsid = models.CharField(max_length=32)
    
    