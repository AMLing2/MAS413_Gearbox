import camelot
import pandas as pd
print('Camelot version:', camelot.__version__)

# Extract tables from the PDF: Single Row Deep Groove Ball Bearings
# start = 262
# end = 309
start = 1
end = 5
tables = camelot.read_pdf('SKFbearings_20250321_ballBearings-1-5.pdf', pages=f'{start}-{end}', flavor='stream')
print('\nBall Bearing Tables found: ', len(tables))

for i, table in enumerate(tables):
    print(f'Saved table {i + start}')
    # Combine all tables into one DataFrame
    combinedDataframe = pd.concat([table.df for table in tables])
    
# Export the combined DataFrame to a single CSV file
combinedDataframe.to_csv('combined_ballBearings.csv')
print('Saved all tables to combined_ballBearings.csv')

# Notes:
    # Prøver å hente ut tabeller for diameterne og fatigue load limit fra SKF katalogen,
    # da vi har tilgang på X og Y for dem i den master dokumentet jeg la på disc.
    # Issuet er at man må sette opp en del greier for å få eksportert riktig, og vi 
    # skal bare ha de tabellene fra annenhver side... Så enten her eller i matlab (test 
    # script i denne folderen) at man fjerner misc, evt endrer settings for Camelot
    # eller bytter til en annen library.. Bruker ghostscript versjon 9.55.0, som må
    # lastes ned lokalt, og brukte conda for å importere Camelot. Om vi gjør det slik,
    # kan det fortsatt automatiseres bearing valg ved å beregne equivalent dynamic load
    # for alle og iterere gjennom bearings. Ellers må vi ta det manuelt.. Går kanskje 
    # det og, muligens sunken cost fallacy @ this point... Tips: test med 1-5 pdf'en først
    # den er mindre så d ikke tar så lang tid.