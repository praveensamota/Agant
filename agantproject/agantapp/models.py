from django.db import models

# Create your models here.
class agant(models.Model):
    pdbid = models.CharField(max_length=200)
    l_id = models.CharField(max_length=100)
    _type = models.CharField(max_length=100)
    uniprotid = models.URLField(max_length=100)
    
    def __str__(self):
        return self.pdbid
    