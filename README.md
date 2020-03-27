# covid19

ChemicalChecker based candidates for COVID19.

## Quickstart

```shell
git clone http://gitlabsbnb.irbbarcelona.org/project-specific-repositories/covid19.git
cd covid19
python scripts/fetch_literature.py
export FLASK_APP=web/app.py
flask run 
```

## Tables provided:

### Literature table
1. mol inchikey (link to CC)
2. mol name
3. evidence score
4. description (wrapped)
5. aggregated reference link

### Candidates table
1. mol inchikey
2. mol name
3. mol categories (multiple choice)
4. score (nr similar)
5. score (nr similar different MoA)
6. more scores

## Remote debug

If you are running the flask debug server on a remote machine you'll need to tunnel the port to you machine.
For this it is convenient to have ssh public key access to the machine where you are running the flask server.


Run on your local machine:
```shell
ssh-keygen
ssh-copy-id -i ~/.ssh/id_rsa user@flaskhost
```

Now you can tunnel the remote localhost to your localhost (and see the page in your browser):
```shell
ssh -L 5000:127.0.0.1:5000 user@flaskhost -fN
```