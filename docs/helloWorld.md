
# Markdown
# hello world
###### hello world
[hello world](https://en.wikipedia.org/wiki/%22Hello,_World!%22_program)


```python
#comments which are not run are denoted with #
#hello world
#try running this cell clicking inside it and pressing shift+enter
#python3
print("hello world")
```

## Quick Jupyter Tips
### You can:
* select how the cell is interpreted in the toolbar at the top, this one says "Markdown"
* view a notebook as a powerpoint style presentation by clicking the "play" button at the top right, with button next to it allowing for editing
* edit text at multiple locations by holding down the "command" key
* highlight and edit blocks of text using the "alt" key
* download the notebook in a variety of ways, including html and presentation html (which preserve interactive elements) and pdf which is static


```python
%lsmagic
#shows you all available magics, these allow you to run many different 
#commands on a single line, such as "load_ext" to load functions or "time" to time your code
#cell magics can interpret many other languages
```




    Available line magics:
    %alias  %alias_magic  %autocall  %automagic  %autosave  %bookmark  %cat  %cd  %clear  %colors  %config  %connect_info  %cp  %debug  %dhist  %dirs  %doctest_mode  %ed  %edit  %env  %gui  %hist  %history  %killbgscripts  %ldir  %less  %lf  %lk  %ll  %load  %load_ext  %loadpy  %logoff  %logon  %logstart  %logstate  %logstop  %ls  %lsmagic  %lx  %macro  %magic  %man  %matplotlib  %mkdir  %more  %mv  %notebook  %page  %pastebin  %pdb  %pdef  %pdoc  %pfile  %pinfo  %pinfo2  %popd  %pprint  %precision  %profile  %prun  %psearch  %psource  %pushd  %pwd  %pycat  %pylab  %qtconsole  %quickref  %recall  %rehashx  %reload_ext  %rep  %rerun  %reset  %reset_selective  %rm  %rmdir  %run  %save  %sc  %set_env  %store  %sx  %system  %tb  %time  %timeit  %unalias  %unload_ext  %who  %who_ls  %whos  %xdel  %xmode
    
    Available cell magics:
    %%!  %%HTML  %%SVG  %%bash  %%capture  %%debug  %%file  %%html  %%javascript  %%js  %%latex  %%perl  %%prun  %%pypy  %%python  %%python2  %%python3  %%ruby  %%script  %%sh  %%svg  %%sx  %%system  %%time  %%timeit  %%writefile
    
    Automagic is ON, % prefix IS NOT needed for line magics.




```python
%%perl
#can use magic to run whole cell as perl
print "hello world"
```

    Operator or semicolon missing before %lsmagic at - line 3.
    Ambiguous use of % resolved as operator % at - line 3.
    Illegal modulus zero at - line 2.



```python
%%bash
#or make cell bash
echo "hello world"
```


```python
#bash commands inside python kernel
!echo "hello world"
```


```python
#can even store output of bash commands and use in python
string = !echo "hello world"
print(string)
```

    ['hello world']


## Switching kernels
Jupyter (***Ju***lia***Pyt***hon***R***) is essentially a wrapper of kernels. A kernel takes text, compiles it according their respective rules and outputs machine code which is then executed on the CPU/GPU. Jupyter allows you to switch between different kernels by goint to Kernel-->Change Kernel at the top.


```python
#switch to Bash kernel
echo "hello world"
```


```python
#switch to R kernel
print("hello world")
```


```python
#switch to C kernel
printf("hello world");
```


```python
#switch to python2 kernel
print "hello world"
```

    hello world


# Images
![hello world](https://babyshark.net/wp-content/uploads/2017/03/hello_world1.gif)
<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/2/28/HelloWorld.svg/1200px-HelloWorld.svg.png" style="width: 500px; height: 500px">


```python
#Youtube videos
from datetime import timedelta
from IPython.display import YouTubeVideo
start=int(timedelta(hours=0, minutes=0, seconds=0).total_seconds())
YouTubeVideo("rOU4YiuaxAM", start=start, autoplay=1, theme="light", color="red")
```





        <iframe
            width="400"
            height="300"
            src="https://www.youtube.com/embed/rOU4YiuaxAM?color=red&start=0&theme=light&autoplay=1"
            frameborder="0"
            allowfullscreen
        ></iframe>
        


