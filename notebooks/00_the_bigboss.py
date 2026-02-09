{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "8de13176",
   "metadata": {},
   "outputs": [],
   "source": [
    "import pandas as pd\n",
    "from pathlib import Path\n",
    "\n",
    "meta = pd.read_csv(\"../HumanBreast10X-main/Tables/SampleStats.txt\", sep=\"\\t\")\n",
    "\n",
    "#Match available matrix/barcode files to metadata GEO_IDs\n",
    "data_folder = Path(\"../data/GSE161529_RAW\")\n",
    "\n",
    "rows = []\n",
    "for _, row in meta.iterrows():\n",
    "    geo = row['GEO_ID']\n",
    "    # Find all matching files for this sample\n",
    "    mtx_files = list(data_folder.glob(f\"{geo}_*-matrix.mtx.gz\"))\n",
    "    bc_files = list(data_folder.glob(f\"{geo}_*-barcodes.tsv.gz\"))\n",
    "    for mtx, bc in zip(mtx_files, bc_files):\n",
    "        rows.append({\n",
    "            'GEO_ID': geo,\n",
    "            'MatrixFile': mtx.name,\n",
    "            'BarcodesFile': bc.name,\n",
    "            'SampleName': row['SampleName'],\n",
    "            'Title': row['Title'],\n",
    "            'CellNumAfter': row['CellNumAfter'],\n",
    "            'GenesDetected': row['GenesDetected'],\n",
    "        })\n",
    "\n",
    "summary_table = pd.DataFrame(rows)\n",
    "summary_table.to_csv(\"../results/TheBigBoss.csv\", index=False)\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "03a77c15",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "✓ BigBoss table updated and saved as BigBoss_updated.csv\n"
     ]
    }
   ],
   "source": [
    "#i want to create a paths to load each expression matrxi and barcode programmatically\n",
    "import pandas as pd\n",
    "from pathlib import Path\n",
    "\n",
    "data_folder = Path(\"../data/GSE161529_RAW\")\n",
    "bigboss = pd.read_csv(\"../results/TheBigBoss.csv\")\n",
    "\n",
    "# Helper to find filename by GEO_ID\n",
    "# Assumes each GEO_ID appears only once per sample type\n",
    "\n",
    "def find_file(geo_id, pattern):\n",
    "    files = list(data_folder.glob(f'{geo_id}_*-{pattern}'))\n",
    "    return files[0].name if files else None\n",
    "\n",
    "# Add columns for matrix and barcodes filenames\n",
    "bigboss['MatrixFile'] = bigboss['GEO_ID'].apply(lambda x: find_file(x, 'matrix.mtx.gz'))\n",
    "bigboss['BarcodesFile'] = bigboss['GEO_ID'].apply(lambda x: find_file(x, 'barcodes.tsv.gz'))\n",
    "\n",
    "# Save the updated BigBoss with new columns\n",
    "bigboss.to_csv('../results/TheBigBoss_updated.csv', index=False)\n",
    "print('✓ BigBoss table updated and saved as BigBoss_updated.csv')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "4836550d",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "breast_cancer_scrnaseq",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.10.12"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
