import camelot
print('Camelot version:', camelot.__version__)

# Extract tables from the PDF: Single Row Deep Groove Ball Bearings
start = 262
end = 309
tables = camelot.read_pdf('SKFbearings_20250321.pdf', pages=(start,'-',end), flavor='stream')

print('\nBall Bearing Tables found: ', len(tables))
for i, table in enumerate(tables):
    if (i + 262) % 2 == 0:  # Check if the page number is even
        table.to_csv(f'ballBearing_page_{i + 262}.csv')
        print(f'Saved table from page {i + 262} to CSV')