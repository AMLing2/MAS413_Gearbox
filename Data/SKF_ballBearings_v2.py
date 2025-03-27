import camelot
import pandas as pd
print('Camelot version:', camelot.__version__)

def df_search_text(df,search_text,phrase):
    phrase_found = False
    h_index = None
    for index,line in df.iterrows():
        if line.str.contains(phrase).any():
            phrase_found = True
        if line.str.contains(search_text).any():
            h_index = index
        if phrase_found and h_index is not None:
            break
    if h_index is None:
        print("search text not found")
        exit()
    else:
        return [h_index,phrase_found]

# Extract tables from the PDF: Single Row Deep Groove Ball Bearings
# start = 262
# end = 309
start = 1
end = 47
#tables = camelot.read_pdf('SKFbearings_20250321_ballBearings-1-5.pdf', pages=f'{start}-{end}', flavor='stream')
tables = camelot.read_pdf('SKFbearings_20250321_ballBearings.pdf', pages=f'{start}-{end}', flavor='stream')
print('\nBall Bearing Tables found: ', len(tables))

fulltabs = []
newtab_principle_df = []
newtab_extra_df = []
n = 0
a = 0;
for i, table in enumerate(tables):
        
    print(f'Saved table {i + start}:\n')
    if table.df[0].str.contains("mm").any(): # exclude corrupted tables
        # check table type
        princ_str = "▶"
        [h_index,princ_p] = df_search_text(table.df,"mm",princ_str)
        if princ_p:
            #remove "SKF Explorer bearing" text on some table last rows
            explorer_str = "SKF Explorer bearing"
            if table.df.iloc[-1].str.contains(explorer_str).any():
                newtab = table.df.iloc[:-1]
            else:
                newtab = table.df
            if n == 0:
                newtab_principle_df.append(newtab)
            else: # remove header
                newtab_principle_df.append(newtab.iloc[h_index+1:]) # bad hardcoded
            a += 1

        else: # extra data table
            if n == 0:
                num_padded_rows = newtab_principle_df[0].shape[0] - table.df.shape[0]
                padded_df = pd.DataFrame(
                    # Fill with empty strings for each column:
                    [[''] * len(table.df.columns)] * num_padded_rows,
                    columns=table.df.columns  # Ensure columns match the original DataFrame
                )
                padded_table = pd.concat([padded_df,table.df], ignore_index=True)
                newtab_extra_df.append(padded_table)
            else: # remove header
                newtab_extra_df.append(table.df.iloc[h_index+1:])
            n += 1
            a += 1

# Combine all tables into one DataFrame
principle_combined = pd.concat([df for df in newtab_principle_df])
extra_combined = pd.concat([df for df  in newtab_extra_df])
# fix indexes after concating
principle_combined.reset_index(drop=True, inplace=True)
extra_combined.reset_index(drop=True, inplace=True)
#print(principle_combined.to_string())
#print(extra_combined.to_string())
# combine to a single table:
combinedDataframe = pd.concat([principle_combined,extra_combined], axis=1, ignore_index=False)
    
# Export the combined DataFrame to a single CSV file
combinedDataframe.to_csv('combined_ballBearings.csv',sep=';')
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
