# Chemical Checker based candidates for COVID-19

This is the repository for scripts and web front-end for the COVID-19 project.
The goal is to present a list of bioactive chemical compounds with potential to be effective against COVID-19.

We are actively collecting suggested drugs from the current COVID-19 literature, with different levels of supporting evidence. We then use the Chemical Checker to identify small-molecules with similar chemical and bioactivity features to the reported drugs in a universe of 800 thousand bioactive compounds.


## Tables provided:

### Literature table

A cleaned and tractable version of [literature candidates}(https://docs.google.com/spreadsheets/d/1BesqNdhHoVyldk372JOyNbhKAl0ncVZwloTpWb6YWZw/edit#gid=424137782) where molecules standardized, aggregated and duplications resolved.

It contains the following fields:

1. Molecule InChIKey
2. Molecule name
3. Evidence level
4. Evidence type
5. Mode of action

### Candidates table

1. Molecule InChIKey
2. Molecule name
3. Is Drug?
4. Support score
5. 1e-5 score
6. 1e-4 score
7. 1e-3 score
8. Most similar literature compound
9. Second most similar literature compound
10. Third most similar literature compound


## Quickstart

```shell
git clone http://gitlabsbnb.irbbarcelona.org/project-specific-repositories/covid19.git
cd covid19
python scripts/fetch_literature.py
export FLASK_APP=web/app.py
flask run 
```

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