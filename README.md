# 🧠 CWAS for resting-state fMRI

# Edit run_cwas
1. Change your list of features with the feature name in your halfpipe output
2. Edit all your path with your user ID
3. Change your Case and Control variables (based on your participants.tsv file)
4. Indicate a sequence_col or machine_col if available. Otherwise do not specify these arguments. 

# participants.tsv file
- Make sure to have a participant.tsv file with only GOOD subjects (i.e. subjects who have a good QC)
- Your participants.tsv file need to have the following columns:

| participant_id | diagnosis | gender | age |
| ---------------|-----------|--------|-----|

# Environment
The environment will be automatically created and deleted with the job.
