import csv
from django.core.management.base import BaseCommand
from agantapp.models import agant

class Command(BaseCommand):
    help = 'Import data from CSV file'

    def handle(self, *args, **options):
        file_path = '/home/praveen/agant/agantproject/agantapp/management/commands/output.csv'  # Update this with your actual file path

        with open(file_path, 'r') as file:
            csv_reader = csv.reader(file)
            next(csv_reader)  # Skip header row

            for row in csv_reader:
                pdbid, l_id, _type, uniprotid = row
                agant.objects.create(pdbid=pdbid, l_id=l_id, _type=_type, uniprotid=uniprotid)

        self.stdout.write(self.style.SUCCESS('Data imported successfully'))
