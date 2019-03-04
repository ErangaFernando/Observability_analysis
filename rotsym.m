function R=rot(a,b)
if (a=='x')
    c=cos(b);
    s=sin(b);
    R=[1  0  0;
       0  c -s;
       0  s  c];
end;
if (a=='y')
    c=cos(b);
    s=sin(b);
    R=[ c  0  s;
        0  1  0;
       -s  0  c];
end;
if (a=='z')
    c=cos(b);
    s=sin(b);
    R=[c -s  0;
       s  c  0;
       0  0  1];
end;

   