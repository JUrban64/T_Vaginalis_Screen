# AlphaFold Database - Zpráva o stahování struktur

- **Datum spuštění:** `2026-09-04 12:03:20`
- **Vstupní soubor:** `uniprot_methyltransferases_deduplicated.xlsx`
- **Zvolený formát:** `PDB`
- **Celkem zpracováno proteinů:** **67**

## Souhrnná statistika

| Metrika | Počet | Podíl (%) |
| :--- | :---: | :---: |
| **Úspěšně dostupné / stažené** | **65** | **97.0 %** |
| **Chybějící v AlphaFold DB** | **2** | **3.0 %** |
| **Celkem v seznamu** | **67** | **100.0 %** |

## Distribuce kvality modelů (pLDDT)

- **Průměrné globální pLDDT skóre:** **`88.10`**

| Interval spolehlivosti | Kategorie | Počet struktur | Podíl ze stažených |
| :--- | :--- | :---: | :---: |
| **pLDDT >= 90** | Velmi vysoká spolehlivost (Very high) | 32 | 49.2 % |
| **70 <= pLDDT < 90** | Spolehlivá páteř (Confident) | 32 | 49.2 % |
| **50 <= pLDDT < 70** | Nízká spolehlivost (Low) | 1 | 1.5 % |
| **pLDDT < 50** | Velmi nízká / nestrukturovaná (Disordered) | 0 | 0.0 % |

## Seznam chybějících struktur (nejsou v AlphaFold DB)

Následující proteiny nemají v AlphaFold DB k dispozici 3D model:

| UniProt Entry | Entry Name | Gen | Důvod / Stav | Název proteinu |
| :--- | :--- | :--- | :--- | :--- |
| [`A2EN39`](https://www.uniprot.org/uniprotkb/A2EN39) | A2EN39_TRIV3 | TVAG_087250 | Struktura není v AlphaFold DB dostupná (404 Not Found) | Beige/BEACH domain containing protein |
| [`A2ERV9`](https://www.uniprot.org/uniprotkb/A2ERV9) | A2ERV9_TRIV3 | TVAG_473130 | Struktura není v AlphaFold DB dostupná (404 Not Found) | DUF2428 domain-containing protein |

## Vygenerované soubory

- **Adresář se strukturami:** [`reports_deduplicated/../structures/`](file:///Users/jachymurban/Desktop/TVaginalis_screen/reports_deduplicated/../structures)
- **Chybějící struktury (TSV):** [`reports_deduplicated/missing_structures.tsv`](file://reports_deduplicated/missing_structures.tsv)
- **Chybějící struktury (CSV):** [`reports_deduplicated/missing_structures.csv`](file://reports_deduplicated/missing_structures.csv)
- **Kompletní přehled (TSV):** [`reports_deduplicated/download_summary.tsv`](file://reports_deduplicated/download_summary.tsv)
- **Kompletní přehled (CSV):** [`reports_deduplicated/download_summary.csv`](file://reports_deduplicated/download_summary.csv)
