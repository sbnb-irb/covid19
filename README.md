# covid19

ChemicalChecker based candidates for COVID19.

## Running the debug webapp

```shell
export FLASK_APP=web/app.py
flask run 
```

Tables provided:

## Literature table
1. mol inchikey (link to CC)
2. mol name
3. evidence score
4. description (wrapped)
5. aggregated reference link

## Candidates table
1. mol inchikey
2. mol name
3. mol categories (multiple choice)
4. score (nr similar)
5. score (nr similar different MoA)
6. more scores