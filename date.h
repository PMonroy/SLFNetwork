#ifndef DATE
#define DATE

#define SECONDS_DAY 86400.0

typedef struct{
  unsigned int year; 
  unsigned int month, day;
} date;

#define TIME_TO_DATE(date,time)	     \
  date.year = (time) / 360;		     \
  date.month = ((time) % 360) / 30 + 1;    \
  date.day = ((time) % 360) % 30 + 1; 

#define DATE_TO_TIME(date) (date.year*360)+((date.month-1)*30)+(date.day-1);

#endif
