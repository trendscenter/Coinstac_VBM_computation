import contextlib,traceback,re,os,shutil

@contextlib.contextmanager

def make_file_output(write_dir, template_dict, covariates):

    spm12_types = ['Re','c1Re','c2Re','c3Re','c4Re','c5Re','c6Re','mwc1Re','mwc2Re',
    'mwc3Re','mwc4Re','mwc5Re','mwc6Re','smwc1Re','smwc2Re','smwc3Re','smwc4Re',
    'smwc5Re','smwc6Re','swc1Re','swc2Re','swc3Re','swc4Re','swc5Re','swc6Re','wc1Re',
    'wc2Re','wc3Re','wc4Re','wc5Re','wc6Re']

    for type in spm12_types:

        mode = 0o777

        basepath = os.path.join(os.path.dirname(write_dir),"vbm_outputs")

        covpath = os.path.join(basepath,'covariates')
        if os.path.isdir(covpath) == False:
            os.mkdir(covpath, mode)

        path = os.path.join(basepath,'covariates',type)
        if os.path.isdir(path) == False:
            os.mkdir(path, mode)

        fpath = os.path.join(path,"covariates-"+type+".txt")
        if os.path.isfile(fpath) == False:
            file = open(fpath,"a+")
        else:
            continue

        for i, subj in enumerate(covariates, start=0):

          file_ext = re.search(r".[0-9a-z]+$", subj, re.MULTILINE).group()

          if file_ext == '.gz':
              if "/" in subj or "\\" in subj:
                  #If subject file strings have forward or back slashes
                  subject_str = re.sub(r".*[\\\/]{1}([\w]*).{1}[a-z]*.{1}[a-z]*$", r"\1", subj, 0, re.DOTALL).strip()
              else:
                  #Otherwise
                  subject_str = re.sub(r"([\w]*).{1}[a-z]*.{1}[a-z]*", r"\1", subj, 0, re.DOTALL).strip()

          if file_ext == '.nii':
              if "/" in subj or "\\" in subj:
                  #If subject file strings have forward or back slashes
                  subject_str = re.sub(r".*[\\\/]{1}([\w]*).{1}[a-z]*$", r"\1", subj, 0, re.DOTALL).strip()
              else:
                  #Otherwise
                  subject_str = re.sub(r"([\w]*).{1}[a-z]*", r"\1", subj, 0, re.DOTALL).strip()

          #get src file
          src = os.path.join(basepath,subject_str,'anat','vbm_spm12',type+'.nii')
          #Copy src to dst. (cp src dst)
          filestr = subject_str+'-'+type+'.nii'
          newfile = os.path.join(path,subject_str+'-'+type+'.nii')
          shutil.copy(src, newfile)

          values = []
          row = covariates[subj]

          keys = list(row.keys())
          keys.insert(0, "filename")

          header = ', '.join(map(str, keys))

          if i == 0:
              file.write(header+"\r\n")

          for count, key in enumerate(row):
              values.insert(count, row[key])

          values.insert(0, filestr)
          values = ', '.join(map(str, values))
          file.write(values+"\r\n")

        file.close()
