# ~/.bashrc: executed by bash(1) for non-login shells.

# see /usr/share/doc/bash/examples/startup-files (in the package bash-doc)
# for examples

# If not running interactively, don't do anything
case $- in
    *i*) ;;
      *) return;;
esac

# don't put duplicate lines or lines starting with space in the history.
# See bash(1) for more options
HISTCONTROL=ignoreboth

# append to the history file, don't overwrite it
shopt -s histappend

# for setting history length see HISTSIZE and HISTFILESIZE in bash(1)
HISTSIZE=100000
HISTFILESIZE=20000
hopt -s histappend

# check the window size after each command and, if necessary,
# update the values of LINES and COLUMNS.
shopt -s checkwinsize

# set variable identifying the chroot you work in (used in the prompt below)
if [ -z "${debian_chroot:-}" ] && [ -r /etc/debian_chroot ]; then
    debian_chroot=$(cat /etc/debian_chroot)
fi

# enable color support of ls and also add handy aliases
if [ -x /usr/bin/dircolors ]; then
    test -r ~/.dircolors && eval "$(dircolors -b ~/.dircolors)" || eval "$(dircolors -b)"
    alias ls='ls --color=auto'
    #alias dir='dir --color=auto'
    #alias vdir='vdir --color=auto'

    alias grep='grep --color=auto'
    #alias fgrep='fgrep --color=auto'
    #alias egrep='egrep --color=auto'
fi

# colored GCC warnings and errors
export GCC_COLORS='error=01;31:warning=01;35:note=01;36:caret=01;32:locus=01:quote=01'

# enable programmable completion features (you don't need to enable
# this, if it's already enabled in /etc/bash.bashrc and /etc/profile
# sources /etc/bash.bashrc).
if ! shopt -oq posix; then
  if [ -f /usr/share/bash-completion/bash_completion ]; then
    . /usr/share/bash-completion/bash_completion
  elif [ -f /etc/bash_completion ]; then
    . /etc/bash_completion
  fi
fi

export LC_ALL=""
export LC_CTYPE="en_US.UTF-8"
export TERM=xterm

# Alias
alias lls='ls -lah'
alias ttt='tmux a -t turri'
alias pin2='ssh -i ~/.ssh/id_rsa jmgps@156.35.56.116'
alias labo='ssh -i ~/.ssh/id_rsa labo@156.35.56.78'
alias geo='ssh -i ~/.ssh/id_rsa george@156.35.56.116'
alias deg='ssh jmgps@156.35.56.120'
alias path='readlink -f'
alias sspace_long='perl /home/jmgps/jmgps/downloads/SSPACE-LongRead_v1-1/SSPACE-LongRead.pl'
alias tbl2asn='/data/programs/tbl2asn'

# Export
#expor tPS1="\e[0;35m[\u@\e[m\e[1;35m\h\e[m\e[0;35m:\W]\e[m\$ "
export PS1='\[\033[0;35m\]\u@\[\033[1;35m\]\h\[\033[0;35m\]:\W\[\033[0m\]\$ '
export PATH="/data/programs/mummer-4.0.0beta2/:$PATH"
export PATH="/data/programs/MIX-master/bin/:$PATH"
export PATH="$PATH:/data/programs/sratoolkit.2.9.0-centos_linux64/bin/"
export PATH="$PATH:/data/programs/cufflinks/src/"
export PATH="$PATH:/data/programs/SPAdes-3.12.0-Linux/bin/"
export PATH="$PATH:/data/programs/Quake/bin/"
export PATH="$PATH:/data/programs/trimmomatic/"
export PATH="$PATH:/data/programs/kallisto_linux-v0.44.0/"
export PATH="$PATH:/data/programs/CISA1.3/"
export PATH="$PATH:/data/programs/stringtie-1.3.4d/"
export PATH="$PATH:/usr/local/bin/miniconda43/bin/"
export PATH="$PATH:/mnt/nas1/jmgps/downloads/proovread/bin/"
export PATH="$PATH:/home/jmgps/jmgps/downloads/SSPACE-LongRead_v1-1/"

#A command to see failed disks
file="/home/.disks_comprobations_tmp"
if [ -f $file ]; then
    result=$(grep 'Failed' $file)
    if [ "$result" ]; then
        echo "$(tput setaf 1) $(tput bold) $result $(tput sgr0)" | sed 's/ *phy/--> phy/g'
    fi
fi

# A command to auto-extract files
extract () {
   if [ -f $1 ] ; then
       case $1 in
        *.tar.bz2)      tar xvjf $1 ;;
        *.tar.gz)       tar xvzf $1 ;;
        *.tar.xz)       tar Jxvf $1 ;;
        *.bz2)          bunzip2 $1 ;;
        *.rar)          unrar x $1 ;;
        *.gz)           gunzip $1 ;;
        *.tar)          tar xvf $1 ;;
        *.tbz2)         tar xvjf $1 ;;
        *.tgz)          tar xvzf $1 ;;
        *.zip)          unzip $1 ;;
        *.Z)            uncompress $1 ;;
        *.7z)           7z x $1 ;;
        *)              echo "don't know how to extract '$1'..." ;;
       esac
   else
       echo "'$1' is not a valid file!"
   fi
}
