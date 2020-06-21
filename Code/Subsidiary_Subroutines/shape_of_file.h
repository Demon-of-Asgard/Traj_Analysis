dimf shape_of_file(char *Filename){
  FILE *fr;
  char ch,ch_last = ' ';
  int i,j, nrows = 0,elements_in_file = 0;
  int count;
  double tmp;
  dimf d;

  fr = fopen(Filename, "r");
  if(fr == NULL){
    printf("unable to open %s inside shape_of_file, exiting.\n", Filename);
    exit(0); 
  }
  count = 0;    
  while((fscanf(fr,"%lf",&tmp))!=EOF){
    count ++;
  }
  fclose(fr);

  elements_in_file = count;

  fr = fopen(Filename, "r");
  // counting #rows
  for (ch = fgetc(fr); ch != EOF; ch = fgetc(fr)){
    if ((ch == '\n') && (ch_last != '\n')){
      nrows++;
    }
    ch_last = ch;
  }
  fclose(fr);

  d.rows = nrows;
  d.cols = elements_in_file/nrows;
  //printf("#rows: %d #cols: %d in %s\n",d.rows, d.cols, Filename);
  return(d);
}